# ForgingWeapons
Codes and dataset of the paper entitled "<strong>Fight intensity correlates with stronger and more mechanically efficient weapons in three species of *Aegla* crabs</strong>", published in *Behavioral Ecology and Sociobiology*. 

If you use anything from this repository, please cite the paper:

***

Palaoro, AV; Peixoto, PEC; Benso-Lopes, F; Boligon, DS & Santos, S (2020) Fight intensity correlates with stronger and more mechanically efficient weapons in three species of Aegla crabs. Behavioral Ecology and Sociobiology. DOI: 10.1007/s00265-020-02834-z


-----------------------
#### File summary:</br>

<b>Biomech-GMM-Analyses</b>: *R code*. Code that organizes the data in the .xlsx file, loads the TPS files and run all analyses contained in the paper</br>

<b>males-morpho-data.xlsx</b>: *dataset*. Dataset with all the data used in the analyses of the manuscript, except fofr the shape files.</br>

<b>male-aeglids.tps</b>: *shape file*. File containing the 17 landmarks (6 landmarks and 11 semi-landmarks) digitized on the claws. The 
order of the individuals in this file is the same that you will find in the dataset. For more information on the location of the landmarks and semi-landmarks. Please see Figure 1 in the paper.  


###### Dataset metadata: </br>


| Measures     | Unit        | Legend                                                                 |
| ------------ | ----------- | ---------------------------------------------------------------------- |
| ceph.length  | mm          | Cephalothorax length                                                   |
| claw.length  | mm          | Propodus length                                                        |
| claw.height  | mm          | Propodus height                                                        |
| in.lever     | mm          | Distance from the fulcrum to the base of the dactyl                    |
| out.lever1   | mm          | Distance from the fulcrum to the first tubercle on the dactyl          |
| out.lever2   | mm          | Distance from the fulcrum to the tip of the dactyl                     |
| ma1          | -           | Mechanical advantage = in.lever/out.lever1                             |
| ma2          | -           | Mechanical advantage = in.lever/out.lever2                             |
| maT          | -           | Maximum mechanical advantage                                           |
| apodeme      | mm-squared  | Area of the cuticle in which the muscle attachs to close the claw      |
| icf          | -           | Index of closing force = Apodeme times the mechanical advantage        |
| ID           | -           | Number of the individual                                               |
| name         | -           | Code on the tag of the individual in the photo. Used in the TPS file   |
| sex          | -           | Only males                                                             |
| species      | -           | Three levels. Either *Aegla abtao*, *Aegla longirostri*, or *Aegla denticulata* |
