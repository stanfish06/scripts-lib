#/bin/bash
uv --project ~/Git/bioinfoLib/ run pyGenomeTracks --tracks pyGenomeTracksConfig.ini --region chr8:21,896,000-21,912,000 --trackLabelFraction 0.2 --fontSize 20  --dpi 300 -o fgf17-mesendo.png
uv --project ~/Git/bioinfoLib/ run pyGenomeTracks --tracks pyGenomeTracksConfig.ini --region chr10:103,527,866-103,545,000 --trackLabelFraction 0.2 --fontSize 20  --dpi 300 -o fgf8-mesendo.png
uv --project ~/Git/bioinfoLib/ run pyGenomeTracks --tracks pyGenomeTracksConfig.ini --region chr11:69,580,000-69,596,000 --trackLabelFraction 0.2 --fontSize 20  --dpi 300 -o fgf4-mesendo.png
uv --project ~/Git/bioinfoLib/ run pyGenomeTracks --tracks pyGenomeTracksConfig.ini --region chr17:44,837,000-44,980,000 --trackLabelFraction 0.2 --fontSize 20  --dpi 300 -o wnt3-mesendo.png

uv --project ~/Git/bioinfoLib/ run pyGenomeTracks --tracks pyGenomeTracksConfig_h64.ini --region chr8:21,896,000-21,912,000 --trackLabelFraction 0.2 --fontSize 20  --dpi 300 -o fgf17-pluri.png
uv --project ~/Git/bioinfoLib/ run pyGenomeTracks --tracks pyGenomeTracksConfig_h64.ini --region chr10:103,527,866-103,545,000 --trackLabelFraction 0.2 --fontSize 20  --dpi 300 -o fgf8-pluri.png
uv --project ~/Git/bioinfoLib/ run pyGenomeTracks --tracks pyGenomeTracksConfig_h64.ini --region chr11:69,580,000-69,596,000 --trackLabelFraction 0.2 --fontSize 20  --dpi 300 -o fgf4-pluri.png
uv --project ~/Git/bioinfoLib/ run pyGenomeTracks --tracks pyGenomeTracksConfig_h64.ini --region chr17:44,837,000-44,980,000 --trackLabelFraction 0.2 --fontSize 20  --dpi 300 -o wnt3-pluri.png