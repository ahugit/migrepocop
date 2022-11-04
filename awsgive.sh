#!/bin/sh

keydir=/users/ahu/Dropbox/mykei
sourcefiles=/users/ahu/Dropbox/migrepo
data1=/users/ahu/Dropbox/migrepo/familymigpsid.txt
data2=/users/ahu/Dropbox/data2022/familymig2022.txt
data3=/users/ahu/Dropbox/data2022/taxsim/taxes1983.txt

echo "instance info:"
read instanceid
echo "SEND FILES TO $instanceid"
scp -i $keydir/"key030322.pem" $sourcefiles/*.f90 ec2-user@$instanceid:/home/ec2-user/work/  
scp -i $keydir/"key030322.pem" $sourcefiles/bp093022.txt ec2-user@$instanceid:/home/ec2-user/work/  
scp -i $keydir/"key030322.pem" $sourcefiles/compaws ec2-user@$instanceid:/home/ec2-user/work/  
scp -i $keydir/"key030322.pem" $data1 ec2-user@$instanceid:/home/ec2-user/work/  
scp -i $keydir/"key030322.pem" $data2 ec2-user@$instanceid:/home/ec2-user/work/  
scp -i $keydir/"key030322.pem" $data3 ec2-user@$instanceid:/home/ec2-user/work/  
#read remark
#echo "I am $remark too!"

#scp -i "key030322.pem" $sourcefiles/*.f90 ec2-user@ec2-3-5-223-160.compute-1.amazonaws.com:/home/ec2-user/work/  
