import pandas as pd 
df=pd.read_csv('/home/ubuntu/files/trainingfiles/20_trainingblock1.vcf', header=95, sep='\t')
dropcol=df.columns[0:9]
del df['#CHROM']
dfpos=df['POS'] #so we can add this back in later
df=df.drop(dropcol, axis=1)
df.to_csv("~/testpy/rupasandbox/phasing__file", index=False, sep="\t", line_terminator="\n")
cat phasing__file | \
    awk -F$'\t' '{
        if (NR==1){
            print $0
        }
        else if (NR>1) {
            for (i=2; i<NF; i++){printf substr($i,0,1)","substr($i,3,1)"\t"};
            { printf substr($NF,0,1)","substr($NF,3,1)}; 
            { printf "\n" } 
            }
    }'> phasing__file_pairs.csv
df_arrays=np.array(df)

correct_df_array=df_arrays.transpose()
correct_df_array.shape

"""
need to enter in each pair as a set instead of as a string of o or 1

"""
newarrayshape=correct_df.shape
newdf=np.zeros(list(newarrayshape) + [2])
for i, element in enumerate(correct_df):
    for j, w in enumerate(element):
        newdf[i][j]=w.split(',')
