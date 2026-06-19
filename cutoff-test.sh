#!/bin/bash
gmx make_ndx -f rna_frame10.gro -o index-rna.ndx << EOF
a P OP1 OP2 O5' C5' C4' O4' C3' O3' O2' C2' C1'
name 7 RNA_backbone
EOF
gmx trjconv -f frame10.gro -s frame10.gro -o rna_frame10.gro -n index.ndx << EOF
rna_only
EOF
gmx trjconv -s frame10.gro -f high-temp.xtc -o rna_no_solvent.xtc -pbc mol -center -n index.ndx -ur compact << EOF
rna_only
rna_only
EOF
# Test cutoffs from 2.0 nm to 6.25 nm
for cutoff in 0.1 0.2 0.4 0.6 0.8 1.0 2.0 3.0 3.5 4.0 4.25 4.5 5.0 5.25 5.5 6.0 6.25 ; do
	echo -e "RNA_backbone\nRNA_backbone" | gmx cluster \
		-n index-rna.ndx \
		-cutoff ${cutoff} \
		-f rna_no_solvent.xtc \
		-s rna_frame10.gro \
		-method gromos \
		-o rmsd-clusters-${cutoff}nm.xpm \
		-g cluster-${cutoff}nm.log \
		-dist rmsd-matrix-${cutoff}nm.xvg \
		-ev eigen-values-${cutoff}nm.xvg \
		-sz cluster-sizes-${cutoff}nm.xvg \
		-tr transition-matrix-${cutoff}nm.xpm \
		-ntr transitions-${cutoff}nm.xvg \
		-clid cluster-id-${cutoff}nm.xvg \
		-cl central-structures-${cutoff}nm.pdb
	done

	#============================================================
	# Print the table header cleanly
        #============================================================

	echo "================================================================"
	echo -e "             GROMACS CLUSTERING SUMMARY TABLE"
	echo "================================================================"
	printf "%-15s %-15s %-18s %-15s\n" "Cutoff (nm)" "Cutoff (Å)" "Total Clusters" "Cluster 1 Size"
	echo "----------------------------------------------------------------"

	# Loop through all matching logs sorted numerically by the cutoff value
	for logfile in $(ls -v cluster-*nm.log 2>/dev/null); do
		# 1. Extract cutoff value from the filename using sed
		cutoff_nm=$(echo "$logfile" | sed -E 's/cluster-([0-9.]+)nm\.log/\1/')
		# 2. Convert nanometers to Angstroms using awk
		cutoff_ang=$(awk "BEGIN {print ${cutoff_nm} * 10}")
		# 3. Extract the total number of clusters found
		total_clusters=$(grep "Found" "$logfile" | grep "clusters" | awk '{print $2}')
		# 4. Extract total trajectory frames analyzed
		total_frames=$(grep "Number of structures for matrix" "$logfile" | awk '{print $6}')
		# 5. Extract the size of Cluster 1 (the first cluster line)
		# This finds the row starting with '  1 |' and gets the first number after the pipe
		c1_size=$(grep -E '^[[:space:]]+1[[:space:]]+\|' "$logfile" | head -n 1 | awk '{print $3}')
		# 6. Calculate Cluster 1 Percentage if values exist
		if [ -n "$c1_size" ] && [ -n "$total_frames" ]; then
			c1_pct=$(awk "BEGIN {printf \"%.1f%%\", (${c1_size} / ${total_frames}) * 100}")
			else
				c1_pct="N/A"
				fi
				# Print the row aligned neatly to match the columns
				printf "%-15s %-15s %-18s %-15s\n" "$cutoff_nm" "$cutoff_ang" "$total_clusters" "$c1_pct"
			done
			echo "================================================================"
