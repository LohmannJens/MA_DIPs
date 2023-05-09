#!/bin/bash

for d in Alnaji2019_Cal07 Alnaji2019_NC Alnaji2019_BLEE Alnaji2019_Perth Kupke Pelz Alnaji2021
do
    echo $d
    echo `date +"%Y-%m-%d %T"`
    python classifier.py -n 2 -d $d -g
    mv ../../results/ML/2023-* ../../results/ML/$d
done

echo `date +"%Y-%m-%d %T"`
echo Alnaji2021 n=3
python classifier.py -n 3 -d Alnaji2021 -g
mv ../../results/ML/2023-* ../../results/ML/Alnaji2021_n3
echo `date +"%Y-%m-%d %T"`

echo `date +"%Y-%m-%d %T"`
echo Alnaji2019_IAV
python classifier.py -n 2 -d Alnaji2019_NC,Alnaji2019_Perth,Alnaji2019_Cal07 -g -r
mv ../../results/ML/2023-* ../../results/ML/Alnaji2019_IAV
echo `date +"%Y-%m-%d %T"`

echo `date +"%Y-%m-%d %T"`
echo Alnaji2021 Pelz
python classifier.py -n 2 -d Alnaji2021 -v Pelz -g -r
mv ../../results/ML/2023-* ../../results/ML/Alnaji2021_Pelz
echo `date +"%Y-%m-%d %T"`

echo `date +"%Y-%m-%d %T"`
echo Alnaji2021 Features
python classifier.py -n 2 -d Alnaji2021 -m "feature_comparision"
mv ../../results/ML/2023-* ../../results/ML/Alnaji2021_features
echo `date +"%Y-%m-%d %T"`

echo `date +"%Y-%m-%d %T"`
echo Alnaji2021 Alnaji2019Pelz
python classifier.py -n 2 -d Alnaji2021 -v Alnaji2019_NC,Alnaji2019_Perth,Alnaji2019_Cal07 -g -r
mv ../../results/ML/2023-* ../../results/ML/Alnaji2021_Alnaji2019_IAV
echo `date +"%Y-%m-%d %T"`

echo `date +"%Y-%m-%d %T"`
echo Alnaji2021 shap values
python classifier.py -n 2 -d Alnaji2021 -m shap
mv ../../results/ML/2023-* ../../results/ML/Alnaji2021_shap
echo `date +"%Y-%m-%d %T"`

echo `date +"%Y-%m-%d %T"`
echo Alnaji2021 shap values
python classifier.py -n 2 -d Alnaji2021 -m rule_based_validation
mv ../../results/ML/2023-* ../../results/ML/Alnaji2021_rule_based_validation
echo `date +"%Y-%m-%d %T"`

