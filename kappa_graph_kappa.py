from bioruler.library.kappa_translator import KappaExporter, KappaImporter

import argparse
parser = argparse.ArgumentParser(description='Test translator')
parser.add_argument('-pars', dest='pars', action='store', type=str,
                    help='parser path')
parser.add_argument('-p', dest='p', action='store_const', const=True,
                    default=False, help='print results')
args = parser.parse_args()

in_file = open('test.ka', 'r')

imported = KappaImporter.uncompile(["test.ka"], parser=args.pars)

compiled = KappaExporter.compile_nugget_list(imported)

if args.p:
    f = open('res.txt', 'a+')
    print("---------------------------- Original\n\n%s" % in_file.read(), file=f)

    print("---- Result :\n\n"+"\n\n".join(compiled)+"\n---------------------------- End\n\n", file=f)

else:

    print("---------------------------- Original\n\n%s" % in_file.read())

    print("---- Result :\n\n"+"\n\n".join(compiled)+"\n---------------------------- End\n\n")
