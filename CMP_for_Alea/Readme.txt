Modifications in event file:
	Remove the initial rows in the event file: the ones with code 41

	Delete 'comments' column in the event file because it is causing problems while reading file in matlab

	On one timestamp, there will be only 1 trigger. So, please keep in what follows after what. 
	
	Have a look at the event files and Compassion_TriggerCodes.xlsx. You will understand the modifications

How to do run SOBI:
	Fill the file 'CMP_SOBIfiles.xlsx' to specify on which files you want to run SOBI. Maintain the syntax
	Go to the folder containg matlab scripts, the type 'tSOBI_CMP_wrapper('C:\Alea\CMP_SOBIfiles.xlsx')' in command window
	
How to do SOBI reconstruction
	Fill the file 'CMP_reconSOBIfiles.xlsx' to specify on which files you want to run SOBI reconstruction. Maintain the syntax
	Go to the folder containg matlab scripts, the type 'reconEEGdata_contCMP_wrapper('C:\Alea\CMP_reconSOBIfiles.xlsx')' in command window