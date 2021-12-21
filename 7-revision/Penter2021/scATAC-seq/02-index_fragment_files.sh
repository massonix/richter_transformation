# This script indexes the fragments file, which is necessary for Signac
# It also counts the number of fragments per cell so we can call cells
module load gcc/4.9.3-gold zlib/1.2.8 HTSLIB/latest
dir=/scratch/devel/rmassoni/RICHTER/current/data/Penter2021/scATAC-seq
cd $dir
for i in $(ls); do
  echo $i;
  frag_file=$(ls $i | grep fragments.tsv.gz$);
  zcat "${i}/${frag_file}" | awk -v frag="$i" 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,frag"_"$4,$5}' >| "${i}/${i}_fragments.tsv"
  bgzip "${i}/${i}_fragments.tsv"
  tabix -p bed "${i}/${i}_fragments.tsv.gz";
  zcat "${i}/${i}_fragments.tsv.gz" | cut -f4 | sort | uniq -c | awk '{print $2, $1}' > "${i}/fragments_per_cell.txt"
done

