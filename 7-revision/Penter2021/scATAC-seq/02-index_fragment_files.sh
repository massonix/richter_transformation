# This script indexes the fragments file, which is necessary for Signac
# It also counts the number of fragments per cell so we can call cells
module load gcc/4.9.3-gold zlib/1.2.8 HTSLIB/latest
dir=/scratch/devel/rmassoni/RICHTER/current/data/Penter2021/scATAC-seq
cd $dir
for i in $(ls); do
  echo $i;
  frag_file=$(ls $i | grep fragments.tsv.gz$);
  tabix -p bed "${i}/${frag_file}";
  zcat "${i}/${frag_file}" | cut -f4 | sort | uniq -c | awk '{print $2, $1}' > "${i}/fragments_per_cell.txt"
done

