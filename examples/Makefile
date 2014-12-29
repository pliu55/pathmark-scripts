all : report img

report :
	galaxyPATHMARK.py -z 123454321 -n 10 -b 10 -f "1.0;1.0" data/merge_merged.all.tab data/feature.tab data/pid_500_pathway.tab 

img :
	mkdir img
	circlePlot.py -s data/FOXA1-2_Network.samples -f data/FOXA1-2_Network.features -o 'FOXA1' img data/FOXA1-2_Network.subtype data/merge_merged.all.tab

clean :
	rm -rf signature.tab bootstrap_signature.tab null_signature.tab analysis report img *.log .jobTree*
