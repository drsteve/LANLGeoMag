#This sed line is HIDEOUS, but it fixes a bug in old Doxygen latex templates
#which forces PDF output even when we want DVI
#http://unix.stackexchange.com/questions/1233/replacing-multiple-lines-in-sed-or-awk
1h;1!H;${g;s/\\ifx\\pdfoutput\\undefined\n\\usepackage\[ps2pdf,\n            pagebackref=true,\n            colorlinks=true,\n            linkcolor=blue\n           \]{hyperref}\n\\usepackage{pspicture}\n\\else\n\\usepackage\[pdftex,\n            pagebackref=true,\n            colorlinks=true,\n            linkcolor=blue\n           \]{hyperref}\n\\fi/\\usepackage{ifpdf}\
\\ifpdf\
\\usepackage\[pdftex,\
            pagebackref=true,\
            colorlinks=true,\
            linkcolor=blue,\
            unicode\
           \]{hyperref}\
\\else\
\\usepackage\[ps2pdf,\
            pagebackref=true,\
            colorlinks=true,\
            linkcolor=blue,\
            unicode\
           \]{hyperref}\
\\usepackage{pspicture}\
\\fi/;p;}
