open Yojson.Basic
open Ast

let path = ref ""
let r_name = ref "rules.json"
let a_name = ref "agents.json"
let print = ref false
let kappa_files = ref []

let _ = Arg.parse
    [("-path", Arg.Set_string path, "output path of JSON files");
     ("-a", Arg.Set_string a_name, "name of agents file");
     ("-r", Arg.Set_string r_name, "name of rules file");
     ("-p", Arg.Set print, "print result")]
     (fun x -> kappa_files := x::(!kappa_files))
     "Convert input kappa files into two JSON formated files representing agents and rules"

let eval_link l = match fst l with
  | LNK_VALUE (i,_) -> `Int i
  | LNK_ANY -> `String "?"
  | LNK_SOME -> `String "_"
  | LNK_TYPE (s1, s2) -> `List [`String (fst s1); `String (fst s2)]
  | FREE -> `Null

let rec list_of_internal l = match l with
  | [] -> []
  | (s,l)::q -> (`String s)::(list_of_internal q)

let rec eval_ports l = match l with
  | [] -> []
  | p::q -> (fst p.port_nme, `List (list_of_internal p.port_int))::
            (eval_ports q)

let rec eval_rule_ports l = match l with
  | [] -> []
  | p::q -> (fst p.port_nme, `Assoc ["state", `List (list_of_internal p.port_int);
                                     "binding", eval_link p.port_lnk])::
            (eval_rule_ports q)

let eval_agent a = (fst(fst a), `Assoc (eval_ports (snd a)))

let eval_rule_agent a = `List [`String (fst(fst a)); `Assoc (eval_rule_ports (snd a))]

let rec eval_signatures_rec s = match s with
  | [] -> []
  | a::q -> (eval_agent a)::(eval_signatures_rec q)

let eval_signatures s =  `Assoc (eval_signatures_rec s)

let rec eval_mixture m = match m with
  | [] -> []
  | a::q -> ((eval_rule_agent a)::(eval_mixture q))

let eval_rule r =
  let lhs = `List (eval_mixture r.lhs) in
  let rhs = `List (eval_mixture r.rhs) in

  match r.arrow with
    | RAR -> [`List [lhs; rhs]]
    | LRAR -> [`List [lhs; rhs]; `List [rhs; lhs]]

let rec eval_rules l = match l with
  | [] -> []
  | r::q -> (eval_rule (fst(snd r)))@(eval_rules q)

let () =
  let result : (Ast.agent, Ast.mixture, string, Ast.rule) Ast.compil =
        List.fold_left (KappaLexer.compile Format.std_formatter)
                       Ast.empty_compil (!kappa_files) in

  let agent_decl = eval_signatures (result.signatures) in
  let rules_decl = eval_rules (result.rules) in

  if !path <> "" && not (Sys.file_exists !path) then Unix.mkdir !path 0o777;

  to_file ~len:1024 ~std:true (Filename.concat !path !a_name) agent_decl;
  to_file ~len:1024 ~std:true (Filename.concat !path !r_name) (`List rules_decl);
  if !print then begin
    print_string "Agent declarations :\n";
    print_string (pretty_to_string ~std:true agent_decl);
    print_newline();
    print_string "Rules declarations :\n";
    print_string (pretty_to_string ~std:true (`List rules_decl))
  end
