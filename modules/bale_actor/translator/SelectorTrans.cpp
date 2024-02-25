#include "clang/AST/ASTContext.h"
#include "clang/AST/ASTConsumer.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Frontend/FrontendAction.h"
#include "clang/Rewrite/Core/Rewriter.h"
#include "clang/Tooling/CommonOptionsParser.h"
#include "clang/Tooling/Tooling.h"
#include "llvm/Support/CommandLine.h"

#include <regex>

using namespace clang;
using namespace clang::tooling;
static llvm::cl::OptionCategory MyToolCategory("my-tool options");

using namespace clang::ast_matchers;

#define MAX_DEPTH 2

namespace __internal__ {
    const CXXRecordDecl* searchLambda(const Stmt *s) {
        if (!s) { return nullptr; }
        const LambdaExpr *lambda = dyn_cast<LambdaExpr>(s);
        if (lambda) {
            return lambda->getLambdaClass();
        }
        for (auto const &p : s->children()) {
            auto q =  searchLambda(p);
            if (q) return q;
        }
        return nullptr;
    }

    void findName(const Expr *e, llvm::SmallString<16> &ret) {
        bool done = false;
        const DeclRefExpr *dre = dyn_cast<DeclRefExpr>(e);
        const IntegerLiteral *intL = dyn_cast<IntegerLiteral>(e);
        if (intL) {
            intL->getValue().toStringUnsigned(ret);
            done = true;
        } else {
            const ImplicitCastExpr *ice = dyn_cast<ImplicitCastExpr>(e);
            if (!dre && ice) {
                dre = dyn_cast<DeclRefExpr>(ice->getSubExpr());
                const ImplicitCastExpr *ice2 = dyn_cast<ImplicitCastExpr>(ice->getSubExpr());
                if (ice2) {
                    dre = dyn_cast<DeclRefExpr>(ice2->getSubExpr());
                }
            }
            assert(dre != str::nullptr);
            ret = dre->getDecl()->getName();
            done = true;
        }
        if (!done) {
            llvm::errs() << "Could not get the name of:\n";
            e->dump();
        }
    }

    bool isArgToSend(const VarDecl *d, const Stmt *s) {
        bool ret = false;
        if (!d || !s) {
            return ret;
        }
        const CXXMemberCallExpr *send = dyn_cast<CXXMemberCallExpr>(s);
        if (send) {
            if (send->getMethodDecl()
                && send->getMethodDecl()->getName().compare("send") == 0) {
                const DeclRefExpr *dre = dyn_cast<DeclRefExpr>(send->getArg(1));
                const ImplicitCastExpr *ice = dyn_cast<ImplicitCastExpr>(send->getArg(1));
                if (!dre && ice) {
                    dre = dyn_cast<DeclRefExpr>(ice->getSubExpr());
                }
                if (dre && dre->getDecl() == d) {
                    return true;
                }
            }
        }
        for (auto const &p : s->children()) {
            ret |= isArgToSend(d, p);
        }
        return ret;
    }

    SourceLocation getLocOfNewStatement(ASTContext *Context, const CallExpr *finish, const DeclRefExpr *hs_ptr) {
        const auto& parents = Context->getParents(*finish);
        if (!parents.empty()) {
            const Stmt* parent = parents[0].get<Stmt>();
            if (parent) {
                const auto& gparents = Context->getParents(*parent);
                if (!gparents.empty()) {
                    auto gparent = gparents[0].get<Stmt>();
                    for (auto c: gparent->children()) {
                        const BinaryOperator *bo = dyn_cast<BinaryOperator>(c);
                        if (bo) {
                            const DeclRefExpr *dre = dyn_cast<DeclRefExpr>(bo->getLHS());
                            const CXXNewExpr *cne = dyn_cast<CXXNewExpr>(bo->getRHS());
                            if (dre && cne
                                && dre->getDecl() == hs_ptr->getDecl()) {
                                return dre->getBeginLoc();
                            }
                        }
                    }
                }
            }
        }
    }

    void processCaptures(const CXXRecordDecl *lambda, const unsigned int depth, std::vector<std::vector<const VarDecl *>> &arrays, std::vector<std::vector<const VarDecl *>> &scalars, std::vector<const Stmt*> &lambdas, clang::DiagnosticsEngine &DE) {
        assert(lambda->isLambda() == true && "lambda is expected");

        // Diagnostics
        const unsigned ID_LAMBDA_NONPTR = DE.getCustomDiagID(clang::DiagnosticsEngine::Remark, "[captured] non-pointer var: %0");
        const unsigned ID_LAMBDA_NONPTR_EXCLUDED = DE.getCustomDiagID(clang::DiagnosticsEngine::Remark, "[captured] non-pointer var (excluded because it's an arg to send()): %0");
        const unsigned ID_LAMBDA_PTR = DE.getCustomDiagID(clang::DiagnosticsEngine::Remark, "[captured] pointer var: %0");

        llvm::errs() << "[processCaptures] depth = " << depth << "\n";
        auto const *mDecl = lambda->getLambdaCallOperator();
        lambdas[depth] = mDecl->getBody();
        for (auto const &Capture : lambda->captures()) {
            auto const var = Capture.getCapturedVar();
            if (!var->getType().getTypePtr()->isPointerType()) {
                // if there is send API in lambda avoid args to it
                if(!__internal__::isArgToSend(var, mDecl->getBody())) {
                    DE.Report(lambda->getBeginLoc(), ID_LAMBDA_NONPTR).AddString(var->getName());
                    scalars[depth].push_back(var);
                } else {
                    DE.Report(lambda->getBeginLoc(), ID_LAMBDA_NONPTR_EXCLUDED).AddString(var->getName());
                }
            } else {
                DE.Report(lambda->getBeginLoc(), ID_LAMBDA_PTR).AddString(var->getName());
                arrays[depth].push_back(var);
            }
        }
        auto const *inner_lambda = searchLambda(mDecl->getBody());
        if (inner_lambda) {
            processCaptures(inner_lambda, depth+1, arrays, scalars, lambdas, DE);
        }
#if 0
        if (depth == 0) {
            for (int i = 0; i < MAX_DEPTH; i++) {
                for (auto const &S : scalars[i]) {
                    llvm::errs() << S->getName() << "\n";
                }
            }
        }
#endif
    }

    std::pair<std::string, std::string> synthesizePacketType(const int uniqueID, const int nMBs, std::vector<std::vector<const VarDecl *>> &scalars, std::map<const VarDecl*, std::string> &nameMap, bool &shrank, clang::DiagnosticsEngine &DE) {
        // Diagnostics
        const unsigned ID_PACKET_IS_SCALAR = DE.getCustomDiagID(clang::DiagnosticsEngine::Remark, "packet is scalar");
        const unsigned ID_STRUCT_SHRINKABLE = DE.getCustomDiagID(clang::DiagnosticsEngine::Remark, "packet struct is shrinkable");
        const unsigned ID_STRUCT_NOT_SHRINKABLE = DE.getCustomDiagID(clang::DiagnosticsEngine::Remark, "packet struct is NOT shrinkable");

        if (nMBs == 1 && scalars[0].size() == 1) {
            // single MB & single scalar
            DE.Report(ID_PACKET_IS_SCALAR);
            nameMap.insert(std::pair<const VarDecl*, std::string>(scalars[0][0], "pkt"));
            return std::make_pair(scalars[0][0]->getType().getAsString(), "primitive");
        } else {
            auto packet_info = new std::map<const QualType, unsigned>[nMBs+1];
            std::vector<const VarDecl*> done;
            bool shrinkable = false;
            for (int i = 0; i < nMBs; i++) {
                for (auto const s: scalars[i]) {
                    //
                    QualType type = s->getType();
                    if (packet_info[i].find(type) == packet_info[i].end()) {
                        packet_info[i].insert(std::pair<const QualType, unsigned>(type, 1));
                        packet_info[nMBs].insert(std::pair<const QualType, unsigned>(type, 1));
                    } else {
                        packet_info[i][type]++;
                        if (packet_info[i][type] > packet_info[nMBs][type]) {
                            packet_info[nMBs][type] = packet_info[i][type];
                        }
                    }
                    //
                    if (std::find(done.begin(), done.end(), s) == done.end()) {
                        done.push_back(s);
                    }
                }
            }
            unsigned nSlots = 0;
            for (auto& elem: packet_info[nMBs]) {
                nSlots += elem.second;
            }
            llvm::errs() << "nSlots: " << nSlots << " done.size() = " << done.size() << "\n";
            shrinkable = nSlots < done.size();
            if (shrinkable) {
                shrank = true;
                DE.Report(ID_STRUCT_SHRINKABLE);
                //
                std::string ret = "struct packet" + std::to_string(uniqueID) + "{\n";
                for (auto& elem: packet_info[nMBs]) {
                    const VarDecl*** assign = new const VarDecl**[nMBs];
                    for (int i = 0; i < nMBs; i++) {
                        assign[i] = new const VarDecl*[elem.second];
                        for (int j = 0; j < elem.second; j++) {
                            assign[i][j] = nullptr;
                        }
                    }
                    //
                    llvm::errs() << "[Packet slot assignment]\n";
                    for (int i = 0; i < nMBs; i++) {
                        llvm::errs() << "MB" << std::to_string(i) << "\n";
                        std::vector<const VarDecl*> done;
                        for (auto const s: scalars[i]) {
                            if (elem.first == s->getType()) {
                                for (int j = 0; j < elem.second; j++) {
                                    if ((i - 1 >= 0) && assign[i-1][j] == s) {
                                        assign[i][j] = s;
                                        llvm::errs() << s->getName() << " still in slot " << j << "\n";
                                        done.push_back(s);
                                    }
                                }
                                if (std::find(done.begin(), done.end(), s) == done.end()) {
                                    for (int j = 0; j < elem.second; j++) {
                                        if (assign[i][j] == nullptr) {
                                            assign[i][j] = s;
                                            llvm::errs() << s->getName() << " goes to slot " << j << ", name: " << llvm::Twine("pkt.slot" + std::to_string(j)).str() << "\n";
                                            nameMap.insert(std::pair<const VarDecl*, std::string>(s, llvm::Twine("pkt.slot" + std::to_string(j)).str()));
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // 
                    for (unsigned i = 0; i < elem.second; i++) {
                        ret += elem.first.getAsString() + " " + "slot" + std::to_string(i) + ";\n";
                    }
                }
                ret += "};\n";
                llvm::errs() << ret;
                return std::make_pair("packet" + std::to_string(uniqueID), ret);
            } else {
                DE.Report(ID_STRUCT_NOT_SHRINKABLE);
                std::vector<const VarDecl*> done;
                std::string ret = "struct packet" + std::to_string(uniqueID) + "{\n";
                for (int i = 0; i < nMBs; i++) {
                    for (auto const s: scalars[i]) {
                        if (std::find(done.begin(), done.end(), s) == done.end()) {
                            done.push_back(s);
                            nameMap.insert(std::pair<const VarDecl*, std::string>(s, llvm::Twine("pkt." + s->getName()).str()));
                            ret += s->getType().getAsString();
                            ret += llvm::Twine(" " + s->getName() + ";\n").str();
                        }
                    }
                }
                ret += "};\n";
                llvm::errs() << ret;
                return std::make_pair("packet" + std::to_string(uniqueID), ret);
            }
        }
        return std::make_pair("WRONG", "");
    }

    std::string translateLambdaBody(const Stmt *lambda, std::map<const VarDecl*, std::string> &nameMap, std::string hs_ptr, Rewriter &TheRewriter) {
        assert(lambda->isLambda() == true && "lambda is expected");
        std::string ret = Lexer::getSourceText(CharSourceRange::getTokenRange(lambda->getSourceRange()), TheRewriter.getSourceMgr(), TheRewriter.getLangOpts()).str();

        for (auto const& e: nameMap) {
            // T val -> pkt.val
            std::regex reg_decl(e.first->getType().getAsString() + "\\s+" + e.first->getName().str());
            ret = std::regex_replace(ret, reg_decl, e.second);
#if 1
            std::regex reg_var("([^a-zA-Z0-9_])"+ e.first->getName().str() + "([^a-zA-Z0-9_])");
            ret = std::regex_replace(ret, reg_var, "$1" + e.second + "$2");
#else
            std::regex reg_var2("([(\\[\\s])+"+ e.first->getName().str() + "([(\\]\\s])+");
            ret = std::regex_replace(ret, reg_var2, "$1" + e.second + "$2");
#endif
        }

        // get mailbox id
        std::smatch sm;
        std::regex_search(ret, sm, std::regex("send\\((.+),\\s*(.+),\\s*(.+)\\)"));
        std::string mailbox = sm.str(1);

        // use regex::extended to match multiple lines
        std::regex reg_sel("get_selector.+;", std::regex::extended);
        ret = std::regex_replace(ret, reg_sel, "send(" + mailbox + ", pkt, sender_rank);");

        std::regex reg_sel2(hs_ptr + ".+;", std::regex::extended);
        ret = std::regex_replace(ret, reg_sel2, "send(" + mailbox + ", pkt, sender_rank);");

        //
        std::regex reg_dup("pkt.pkt.");
        ret = std::regex_replace(ret, reg_dup, "pkt.");
        return ret;
    }

}

class SelectorLambdaHandler : public MatchFinder::MatchCallback {
public:
    SelectorLambdaHandler(Rewriter &R) : TheRewriter(R) {};

    virtual void run(const MatchFinder::MatchResult &Result) {
        ASTContext *Context = Result.Context;
        // Diagnostic
        clang::DiagnosticsEngine &DE = Context->getDiagnostics();
        const unsigned ID_LAMBDA_ELIGIBLE = DE.getCustomDiagID(clang::DiagnosticsEngine::Remark, "will be processed (outermost send+lambda)");
        const unsigned ID_LAMBDA_NOT_ELIGIBLE = DE.getCustomDiagID(clang::DiagnosticsEngine::Remark, "will not be processed");

        // hs_ptr
        const auto ActorInstantiation = Result.Nodes.getNodeAs<DeclRefExpr>("hs_ptr");

        // hs_ptr->send();
        const auto ActorSendExpr = Result.Nodes.getNodeAs<CXXMemberCallExpr>("send");

        // hclib::finish
        const auto FinishCallExpr = Result.Nodes.getNodeAs<CallExpr>("finish");

        // selection
        if (!ActorInstantiation || !ActorSendExpr) {
            llvm::errs() << "skipped because ActorInstantiation or ActorSendExpr is NULL\n";
            return;
        }
        auto tname = ActorInstantiation->getDecl()->getType().getAsString();
        if (tname.find("hclib::Actor") == std::string::npos
            && tname.find("hclib::Selector") == std::string::npos) {
            llvm::errs() << "skipped because tname is " << tname << "(not hclib::Actor or hclib::Selector)\n";
            return;
        }

        //
        static std::vector<const Stmt*> processed;

        // lambda
        auto const *Lambda = Result.Nodes.getNodeAs<CXXRecordDecl>("Lambda");
        auto const *mDecl = Lambda->getLambdaCallOperator();

        if (std::find(processed.begin(), processed.end(), mDecl->getBody()) != processed.end()) {
            DE.Report(Lambda->getBeginLoc(), ID_LAMBDA_NOT_ELIGIBLE);
            return;
        } else {
            DE.Report(Lambda->getBeginLoc(), ID_LAMBDA_ELIGIBLE);
        }

        // Analyze the captured vars by the lambda
        // 1. Non pointer type: treated as an index to global arrays
        // 2. Pointer type: treated as a global array
        // TODO: maybe some selection is required
        std::vector<std::vector<const VarDecl *>> arrays(MAX_DEPTH);
        std::vector<std::vector<const VarDecl *>> scalars(MAX_DEPTH);
        std::vector<const Stmt*> lambdas(MAX_DEPTH);

        __internal__::processCaptures(Lambda, 0, arrays, scalars, lambdas, DE);
        processed.insert(processed.end(), lambdas.begin(), lambdas.end());

        // # of mailboxes = depth
        int depth = 0;
        for (; depth < MAX_DEPTH; depth++) {
            llvm::errs() << "[captured level-" << depth << "]\n";
            if (scalars[depth].empty()) {
                break;
            }
            for (auto const &S : scalars[depth]) {
                llvm::errs() << S->getName() << "(type = " << S->getType().getAsString() << ")\n";
            }
        }

        const int nMBs = depth;
        static int uniqueID = 0;
        std::map<const VarDecl*, std::string> nameMap;
        std::pair<std::string, std::string> packet;
        bool shrank = false;
        packet = __internal__::synthesizePacketType(uniqueID, nMBs, scalars, nameMap, shrank, DE);
        std::string packet_type = packet.first;
        std::string packet_def = packet.second;
        llvm::errs() << "packet_type: " << packet_type << "\n";

        // Synthesize a selector class
        ActorInstantiation->getDecl()->dump();
        const VarDecl *ActorInstance = dyn_cast<VarDecl>(ActorInstantiation->getDecl());
        // actorInstance -> uniqueID mapping;
        static std::map<const VarDecl*, int> instanceMap;
        // Synthesize
        if (instanceMap.find(ActorInstance) == instanceMap.end())
        {
            std::string mes = "\n";
            if (packet_def.compare("primitive") != 0) {
                mes += packet_def;
            }
            const std::string CLASS = "SynthesizedActor" + std::to_string(uniqueID);
            mes += "class " + CLASS + " : public hclib::Selector<" + std::to_string(nMBs) + ", " + packet_type + "> { \n";
            mes += "public: \n";
            // output arrays
            // (outermost lambda should captures the all array)
            for (auto const v : arrays[0]) {
                mes += v->getType().getAsString() + " " + v->getName().str() + ";\n";
            }

            // lambda body
            for (unsigned int mb = 0; mb < nMBs; mb++) {
                mes += "void process" + std::to_string(mb) + "(";
                mes += packet_type + " pkt, ";
                mes += "int sender_rank)";
                mes += __internal__::translateLambdaBody(lambdas[mb], nameMap, ActorInstance->getName().str(), TheRewriter);
                mes += "\n";
            }

            // Constructor
            mes += CLASS + "(";
            for (auto const v : arrays[0]) {
                mes += v->getType().getAsString() + " _" + v->getName().str();
                if (v != *(arrays[0].end() - 1)) mes += ", ";
            }
            mes += ")";
            for (auto const v : arrays[0]) {
                if (v == *(arrays[0].begin())) mes += ": ";
                mes += v->getName().str() + "(_" + v->getName().str() + ")";
                if (v != *(arrays[0].end() - 1)) mes += ", ";
                else mes += " ";
            }

            // mb[0].process = [this](pkt_type pkt, int sender_rank) { this->process(pkt, sender_rank);};
            mes += "{\n";
            for (unsigned int mb = 0; mb < nMBs; mb++) {
                mes += "mb[" + std::to_string(mb) + "].process = [this](" + packet_type + " pkt, int sender_rank) { this->process" + std::to_string(mb) + "(pkt, sender_rank); };\n";
            }
            mes += "}\n";
            mes += "};\n";

            // instantiation
            mes += CLASS + " *" + ActorInstance->getName().str();
            mes += " = new " + CLASS + "(";
            for (auto const v : arrays[0]) {
                mes += v->getName().str();
                if (v != *(arrays[0].end() - 1)) mes += ", ";
            }
            mes += ");\n";
            mes += "//";
            //ActorInstantiation->getDecl()->getBeginLoc().print(llvm::errs(), TheRewriter.getSourceMgr());
            if (!ActorInstance->hasGlobalStorage()) {
                TheRewriter.InsertText(ActorInstantiation->getDecl()->getBeginLoc(), mes, true, true);
            } else {
                SourceLocation loc = __internal__::getLocOfNewStatement(Context, FinishCallExpr, ActorInstantiation);
                TheRewriter.InsertText(loc, mes, true, true);
            }
            instanceMap.insert(std::pair<const VarDecl*, int>(ActorInstance, uniqueID));
        }
        // Update actor->send
        {
            // packet
            const std::string PACKET = "pkt" + std::to_string(uniqueID);
            // 
            if (instanceMap.find(ActorInstance) != instanceMap.end()
                && instanceMap[ActorInstance] != uniqueID
                && packet_def.compare("primitive") != 0) {
                packet_type = "packet" + std::to_string(instanceMap[ActorInstance]);
            }
            std::string mes1 = packet_type + " " + PACKET + ";\n";
            if (scalars[0].size() == 1) {
                mes1 += PACKET + " = " + llvm::Twine(scalars[0][0]->getName()).str() + ";\n";
            } else {
                for (auto const &S: scalars[0]) {
                    if (!shrank) {
                        mes1 += llvm::Twine(PACKET + "." + S->getName() + " = " + S->getName() + ";\n").str();
                    } else {
                        std::string tmp = nameMap[S];
                        std::regex reg_pkt("pkt\\.");
                        tmp = std::regex_replace(tmp, reg_pkt, "");
                        mes1 += llvm::Twine(PACKET + "." + tmp + " = " + S->getName() + ";\n").str();
                    }
                }
            }

            // see if the return variable of send() is used
            const VarDecl *sendReturnDecl = nullptr;
            SourceLocation insertPoint;
            const auto& parents = Context->getParents(*ActorSendExpr);
            if (!parents.empty()) {
                const Stmt* parent = parents[0].get<Stmt>();
                if (parent) {
                    const auto gparents = Context->getParents(*parent);
                    if (!gparents.empty()) {
                        const VarDecl* decl  = gparents[0].get<VarDecl>();
                        if (decl) {
                            insertPoint = decl->getBeginLoc();
                            sendReturnDecl = decl;
                        }
                    }
                }
            }
            llvm::SmallString<16> orgArg0;
            llvm::SmallString<16> orgArg1;
            __internal__::findName(ActorSendExpr->getArg(0), orgArg0);
            __internal__::findName(ActorSendExpr->getArg(1), orgArg1);
            std::string mes2 = ActorInstance->getName().str() + "->send(";
            mes2 += orgArg0.str();
            mes2 += ", pkt" + std::to_string(uniqueID) + ", ";
            mes2 += orgArg1.str();
            mes2 += ")";
            if (sendReturnDecl) {
                TheRewriter.ReplaceText(SourceRange(insertPoint, ActorSendExpr->getEndLoc()), llvm::Twine(mes1 + "bool " + sendReturnDecl->getName() + " = " + mes2).str());
            } else {
                TheRewriter.ReplaceText(SourceRange(ActorSendExpr->getBeginLoc(), ActorSendExpr->getEndLoc()), mes1+mes2);
            }
        }
        uniqueID++;
    }

private:
    Rewriter &TheRewriter;
};

class SelectorTransConsumer : public clang::ASTConsumer {
public:
    explicit SelectorTransConsumer(Rewriter &R)
        : HandlerForLambda(R) {

        auto Lambda = expr(hasType(cxxRecordDecl(isLambda()).bind("Lambda")));
        // This matches some internal repl of ->send();
#if 0
        Matcher.addMatcher(cxxMemberCallExpr(on(declRefExpr().bind("hs_ptr")),callee(cxxMethodDecl(hasName("send"))),
                                             hasArgument(2, Lambda)).bind("send"),
                           &HandlerForLambda);
#elif 1
        Matcher.addMatcher(callExpr(callee(functionDecl(hasName("hclib::finish"))),
                                    hasAnyArgument(hasDescendant(
                                                       lambdaExpr(
                                                           forEachDescendant(
                                                               cxxMemberCallExpr(on(declRefExpr().bind("hs_ptr")),callee(cxxMethodDecl(hasName("send"))),
                                                                                 hasArgument(2, Lambda)).bind("send")
                                                               )
                                                           )
                                                       )
                                        )).bind("finish"), &HandlerForLambda);
#else
        Matcher.addMatcher(cxxMemberCallExpr(callee(cxxMethodDecl(hasName("send"))),
                                             hasArgument(1, Lambda),
                                             unless(hasAncestor(cxxMemberCallExpr(callee(cxxMethodDecl(hasName("send"))))))
                               ).bind("send"),
                           &HandlerForLambda);
#endif
    }

    virtual void HandleTranslationUnit(clang::ASTContext &Context) {
        // Run the matchers when we have the whole TU parsed.
        Matcher.matchAST(Context);
    }
private:
    SelectorLambdaHandler HandlerForLambda;
    MatchFinder Matcher;
};

class SelectorTransAction : public clang::ASTFrontendAction {
public:
    void EndSourceFileAction() override {
        SourceManager &SM = TheRewriter.getSourceMgr();
        llvm::errs() << "** EndSourceFileAction for: "
                     << SM.getFileEntryForID(SM.getMainFileID())->getName() << "\n";
        TheRewriter.getEditBuffer(SM.getMainFileID()).write(llvm::outs());
    }

    virtual std::unique_ptr<clang::ASTConsumer> CreateASTConsumer(
        clang::CompilerInstance &Compiler, llvm::StringRef InFile) {
        llvm::errs() << "** Creating AST consumer for: " << InFile << "\n";
        TheRewriter.setSourceMgr(Compiler.getSourceManager(), Compiler.getLangOpts());
        return std::unique_ptr<clang::ASTConsumer>(
            new SelectorTransConsumer(TheRewriter));
    }
private:
    Rewriter TheRewriter;
};

int main(int argc, const char **argv) {
    auto OptionsParserOrError = CommonOptionsParser::create(argc, argv, MyToolCategory);
    if (auto Err = OptionsParserOrError.takeError()) {
        llvm_unreachable("Option Error");
    }
    CommonOptionsParser &OptionsParser = *OptionsParserOrError;
    ClangTool Tool(OptionsParser.getCompilations(),
                   OptionsParser.getSourcePathList());

    return Tool.run(newFrontendActionFactory<SelectorTransAction>().get());
}
