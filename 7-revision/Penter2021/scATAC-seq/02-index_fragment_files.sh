# This script indexes the fragments file, which is necessary for Signac
# It also counts the number of fragments per cell so we can call cells
# module load gcc/4.9.3-gold zlib/1.2.8 HTSLIB/latest
dir=/scratch/devel/rmassoni/RICHTER/current/
for i in $(ls "${dir}/data/Penter2021/scATAC-seq/"); do
  echo $i;
  gse_dir="${dir}/data/Penter2021/scATAC-seq/${i}";
  frag_file="${gse_dir}/$(ls $gse_dir)"
  tabix -p bed $frag_file;
  zcat $frag_file | cut -f4 | sort | uniq -c | awk '{print $2, $1}' > "${gse_dir}/fragments_per_cell.txt"
done

