#!/bin/sh

keydir=/users/ahu/Dropbox/mykei
WHERETO=/users/ahu/Dropbox/output1

#scp -i "key030322.pem" ~/Dropbox/migrepo/*.f90 ec2-user@ec2-3-5-223-160.compute-1.amazonaws.com:/home/ec2-user/awsruns/  
#scp -i "key030322.pem" ~/Dropbox/migrepo/*.f90 ec2-user@ec2-3-5-223-160.compute-1.amazonaws.com:/home/ec2-user/awsruns/  

echo "instance info:"
read instanceid
echo "I will now get files from $instanceid"
scp -i $keydir/"key030322.pem" ec2-user@$instanceid:/home/ec2-user/work/*  $WHERETO/  
echo "Done"


