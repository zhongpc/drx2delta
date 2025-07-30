rm -rf TM_migration_drx;
rm *.cif;

# copy the original label drx structure to the current directory
cd ./original_label_drx;
for file in ./drx_fusion*;
do echo $file;
cp $file ..;
done;
cd ..;

# run the site exchange job
for file in ./drx_fusion*;
do echo $file;
for i in {1..50};  # apply 50 site exchange trials
do echo $i;
python py_a2a_exchange.py -i $file -n $i;
done;
done
