#!/bin/sh

keydir=/users/ahu/Dropbox/mykei
WHERETO=/users/ahu/Dropbox/output2

#scp -i "key030322.pem" ~/Dropbox/migrepo/*.f90 ec2-user@ec2-3-5-223-160.compute-1.amazonaws.com:/home/ec2-user/awsruns/  
#scp -i "key030322.pem" ~/Dropbox/migrepo/*.f90 ec2-user@ec2-3-5-223-160.compute-1.amazonaws.com:/home/ec2-user/awsruns/  
#ssh -i "key030322.pem" ec2-user@ec2-44-212-96-245.compute-1.amazonaws.com

echo "instance info:"
read instanceid
echo "I will now connect to $instanceid"
ssh -i $keydir/"key030322.pem" ec2-user@$instanceid  
echo "Done"


