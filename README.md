# RF4Bioinformatics
Random Forest Wrapper Code for book chapter "random forest for Bioinformatics" 

<BR>

If you use this code, please cite : 

- Major paper: 

@incollection{qi2012random,
  title={Random forest for bioinformatics},
  author={Qi, Yanjun},
  booktitle={Ensemble Machine Learning},
  pages={307--323},
  year={2012},
  publisher={Springer}
}

<BR>

PDF @ http://www.cs.cmu.edu/%7Eqyj/papersA08/11-rfbook.pdf

- Related Papers,  

Y. Qi, HK. Dhiman, et al, Z. Bar-Joseph, J. Klein-Seetharaman,(2009) 
"Systematic prediction of human membrane receptor interactions"
PROTEOMICS 2009, 9, 5243-5255


<BR>

---------------------------------------------

I assume that you have the "g77" command in
your command list. (most unix machines have this installed.)


If not, for windows, you can install the software:
"MinGW".
Remember to add the g77 in your windows command
list after installations. 

---------------------------------------------

It is very easy to run this code: 

1. Just put your parameter in a parameter file
for example: testpara.file
totally 10 parameters (all related files' names 
should be also in the input parameter file.)


2. Then run the perl wrapper
perl change_RF_codeRunParaF.pl testpara.file


---------------------------------------------

I assume that your input feature files have been 
pre-processed ==> which means 
they contain all real features and the features 
have no missing values.


---------------------------------------------

In the subdirectory, there exsit the 
example files configurying in the "testpara.file"

---------------------------------------------

perl wrapper "process_RFtestOut.pl" is an extra script.
Since the RF output contains the voting from 
all trees about postive leaves or negative leaves..
Thus one summary score could be just the 
==> positive vote score - negative vote score

This wrapper would convert the direct output 
RF file into the summary score file as described. 
