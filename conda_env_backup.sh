#! /bin/bash/
backup_date=$(date +%Y%m%d)
new_dir=${backup_date}_conda_env
mkdir new_dir

conda env list|sed -n "4,$ p"|awk '{print$1}'|while read env_name;
do
	#echo $id
	`conda env export --file ${new_dir}/${env_name}_evn_${backup_date}.yml --name ${env_name}`
done
