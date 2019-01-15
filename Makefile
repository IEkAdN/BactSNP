all:
	perl -pe "s|ROOT_PATH_DUMMY|$(PWD)|;" bactsnp.sh >bactsnp
	chmod 755 bactsnp
	make -C src/sam2pseudo;    mv src/sam2pseudo/a                sub_exe/sam2pseudo
	make -C src/delta2pseudo;  mv src/delta2pseudo/a              sub_exe/delta2pseudo
	make -C src/merge_fa;      mv src/merge_fa/a                  sub_exe/merge_fa
	make -C src/fa2snp;        mv src/fa2snp/a                    sub_exe/fa2snp
	make -C src/sam2fq;        mv src/sam2fq/a                    sub_exe/sam2fq
	make -C src/merge_region;  mv src/merge_region/a              sub_exe/merge_region
	make -C src/mask_region;   mv src/mask_region/a               sub_exe/mask_region
	make -C src/snp2fa;        mv src/snp2fa/a                    sub_exe/snp2fa
	make -C src/Platanus_trim; mv src/Platanus_trim/platanus_trim sub_exe/
	make -C src/bwa;           mv src/bwa/bwa                     sub_exe/
	make -C src/Platanus;      mv src/Platanus/platanus           sub_exe/

clean:
	rm bactsnp
	rm sub_exe/sam2pseudo
	rm sub_exe/delta2pseudo
	rm sub_exe/merge_fa
	rm sub_exe/fa2snp
	rm sub_exe/sam2fq
	rm sub_exe/merge_region
	rm sub_exe/mask_region
	rm sub_exe/snp2fa
	rm sub_exe/platanus_trim
	rm sub_exe/bwa
	rm sub_exe/platanus
	rm src/sam2pseudo/*.o
	rm src/sam2pseudo/*.d
	rm src/delta2pseudo/*.o
	rm src/delta2pseudo/*.d
	rm src/merge_fa/*.o
	rm src/merge_fa/*.d
	rm src/fa2snp/*.o
	rm src/fa2snp/*.d
	rm src/sam2fq/*.o
	rm src/sam2fq/*.d
	rm src/merge_region/*.o
	rm src/merge_region/*.d
	rm src/mask_region/*.o
	rm src/mask_region/*.d
	rm src/snp2fa/*.o
	rm src/snp2fa/*.d
	make clean -C src/Platanus_trim
	make clean -C src/bwa
	make clean -C src/Platanus
