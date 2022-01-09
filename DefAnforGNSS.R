################################################################################
#                                                                              #
#     This file contains the code for the  DefAn-for-GNSS-networks.            #
#     DefAn-for-GNSS-networks is created using the RStudio Desktop.            #
#                                                                              #
#     Copyright (C) 2021  Burhaneddin Bilgen                                   #
#                                                                              #
#     This program is free software: you can redistribute it and/or modify     #
#     it under the terms of GNU GENERAL PUBLIC LICENSE. See                    #
#     https://www.gnu.org/licenses/ for terms and conditions.                  #
#                                                                              #
#     This program is distributed in the hope that it will be useful,          #
#     but without any warranty; without even the implied warranty of           #
#     merchantability or fitness for a particular purpose.                     #
#                                                                              #
#                                                                              #
#     Please report any bug, error or suggestion to bbilgen@ktun.edu.tr        #
#                                                                              #
################################################################################



library(dplyr)
library(foreign)
library(stringr)
library(ISwR)
library(matlib)

# Reading the Bernese *.OUT Files
if (interactive() && .Platform$OS.type == "windows"){
  print("Select *.OUT files")
  repeat{
    fils = sort(choose.files(default = "*.OUT", caption = "Select *.OUT files", multi = TRUE), decreasing = FALSE)
    variant = c(paste ("first", seq(1,length(fils)), sep = ""))
    for (h in 1:length(fils)){
      variant [h] = list(variant[[h]][1])
    }
    for (h in 1:length(variant)){
      variant [[h]] = readLines(fils[h])
    }
    if(length(fils) == 1){
      cat(sprintf("At least two *.OUT files have to been selected. You chose %s file.",length(fils)))
      cat(sprintf("\n"))
      rstudioapi::showDialog(title = "Warning!",
                             message = "At least two *.OUT files have to been <b>selected.</b>")
    }
    if(length(fils)>1){
      break
    }
  }
}else {
  repeat{
    data.path = readline(prompt = "Please type the directory of *.OUT files (eg: C:/.../data files/) = ")
    datafiles = list.files(path=data.path, pattern="\\.OUT$")
    fils = sort(paste(data.path, datafiles, sep = ""), decreasing = FALSE)
    variant = c(paste ("first", seq(1,length(fils)), sep = ""))
    for (h in 1:length(fils)){
      variant [h] = list(variant[[h]][1])
    }
    for (h in 1:length(variant)){
      variant [[h]] = readLines(fils[h])
    }
    rm(datafiles, data.path)
    if(length(fils) == 1){
      cat(sprintf("At least two *.SNX files have to been selected. You chose %s file.",length(fils)))
      cat(sprintf("\n"))
    }
    if(length(fils)>1){
      break
    }
  }
}

# Reading the Bernese *.SNX Files
if (interactive() && .Platform$OS.type == "windows"){
  print("Select *.SNX files")
  repeat{
    fils = sort(choose.files(default = "*.SNX", caption = "Select *.SNX files", multi = TRUE), decreasing = FALSE)
    variant2 = c(paste ("firstSNX", seq(1,length(fils)), sep = ""))
    for (h in 1:length(fils)){
      variant2 [h] = list(variant2[[h]][1])
    }
    for (h in 1:length(variant2)){
      variant2 [[h]] = readLines(fils[h])
    }
    if(length(fils) == 1){
      cat(sprintf("At least two *.SNX files have to been selected. You chose %s file.",length(fils)))
      cat(sprintf("\n"))
      rstudioapi::showDialog(title = "Warning!",
                             message = "At least two *.SNX files have to been <b>selected.</b>")
    }
    if(length(fils)>1){
      break
    }
  }
}else {
  repeat{
    data.path = readline(prompt = "Please type the directory of *.SNX files (eg: C:/.../data files/) = ")
    datafiles = list.files(path=data.path, pattern="\\.SNX$")
    fils = sort(paste(data.path, datafiles, sep = ""), decreasing = FALSE)
    variant2 = c(paste ("firstSNX", seq(1,length(fils)), sep = ""))
    for (h in 1:length(fils)){
      variant2 [h] = list(variant2[[h]][1])
    }
    for (h in 1:length(variant2)){
      variant2 [[h]] = readLines(fils[h])
    }
    rm(datafiles, data.path)
    if(length(fils) == 1){
      cat(sprintf("At least two *.SNX files have to been selected. You chose %s file.",length(fils)))
      cat(sprintf("\n"))
    }
    if(length(fils)>1){
      break
    }
  }
}

numextract = function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*[^[:digit:].]")
}

# Taking the degree of freedom, A posteriori RMS of unit weight and cartesian coordinates of stations for all epochs
if (length(variant) == length(variant2)){
  
  #Defining the general parameters and required variables
  options(digits = 12)
  Dof = c(rep(0,length(variant)))
  rms = c(rep(0,length(variant)))
  start1 = c(rep(0,length(variant))) 
  finish1 = c(rep(0,length(variant)))
  x = list()
  pid = list()
  dizi_pid = list()
  b1 = nchar("Estimated value")
  c1 = nchar("RMS error")
  v = list()
  v1 = list()
  epochs = list()
  rnms = list()
  rmse = list()
  rmse1 = list()
  rmse2 = list()
  rnms_rms = list()
  
  # DOF, RMS and cartesian coordinates are taken from whole data in the loop
  for(h in 1:length(variant)){
    Dof [h] = readr::parse_number(variant[[h]][grep("Degree of freedom",variant[[h]])])
    rms [h] = readr::parse_number(variant[[h]][grep("A posteriori RMS of unit weight",variant[[h]])])
    start1 [h] = grep("Reference epoch:",variant[[h]]) + 2
    finish1 [h] = grep("Troposphere parameters:",variant[[h]]) - 4
    x [[h]] = variant [[h]] [start1[h]:finish1[h]]
    start1 [h] = regexpr("Station name", x[[h]][1])
    finish1 [h] = start1[h] + 3
    pid [h]= list(c(substr(x[[h]][3],start1[h],finish1[h])))
    for(i in 4:length(x[[h]])){
      pid [[h]] = rbind (pid [[h]], c(substr(x[[h]] [i], start1 [h], finish1 [h])))
    }
    pid [[h]] = matrix(str_extract(pid [[h]], "[[A-Z]||[:digit:]][[A-Z]||[:digit:]][[A-Z]||[:digit:]][[A-Z]||[:digit:]]"))
    dizi_pid [[h]] = seq(1, length(pid[[h]]), 8)
    pid [[h]] = as.matrix(pid[[h]][dizi_pid[[h]]])
    start1 [h] = regexpr("Estimated value",x[[h]][1])
    finish1 [h] = start1[h]+b1
    v [[h]]= matrix(as.numeric(numextract(substr(x [[h]][3], start1[h], finish1[h]))))
    for(i in 4:length(x[[h]])){
      v [[h]] = rbind(v [[h]], matrix(as.numeric(numextract(substr(x [[h]][i], start1[h], finish1[h])))))
    }
    v [[h]] = na.exclude(v[[h]])
    dizi_x = seq(1,length(v[[h]]),6)
    dizi_y = seq(2,length(v[[h]]),6)
    dizi_z = seq(3,length(v[[h]]),6)
    dizi = c(dizi_x, dizi_y, dizi_z)
    dizi = sort(dizi)
    v1 [[h]] = as.matrix(v [[h]] [dizi])
    epochs [[h]] = v[[h]][dizi_x]
    epochs [[h]] = cbind(epochs[[h]],v[[h]][dizi_y],v[[h]][dizi_z])
    colnames(epochs[[h]]) = c("X","Y","Z")
    rownames(epochs[[h]]) = c(pid[[h]])
    epochs [[h]] = epochs [[h]] [order(row.names(epochs [[h]]), decreasing = FALSE), ]
    rnms [[h]] = row.names(epochs[[h]])
    j = 1
    for (i in seq(1, length(v1 [[h]]), 3)){
      v1 [[h]] [i,1] = epochs [[h]] [rnms [[h]] [j], "X"]
      j = j+1
    }
    j = 1
    for (i in seq(2, length(v1 [[h]]), 3)){
      v1 [[h]] [i,1] = epochs [[h]] [rnms [[h]] [j], "Y"]
      j = j+1
    }
    j = 1
    for (i in seq(3, length(v1 [[h]]), 3)){
      v1 [[h]] [i,1] = epochs [[h]] [rnms [[h]] [j], "Z"]
      j = j+1
    }
    print(epochs[[h]])
    cat(sprintf("\n"))
  }
  
  #RMS errors are taken from whole data in the loop
  for(h in 1:length(variant)){
    start1 [h] = regexpr("RMS error",x[[h]][1])
    finish1 [h] = start1[h]+c1-1
    rmse [[h]] = matrix(as.numeric(substr(x [[h]][3], start1[h], finish1[h])))
    for(i in 4:length(x[[h]])){
      rmse [[h]] = rbind(rmse [[h]], matrix(as.numeric(substr(x [[h]][i], start1[h], finish1[h]))))
    }
    rmse [[h]] = na.exclude(rmse[[h]])
    dizi_x = seq(1,length(rmse[[h]]),6)
    dizi_y = seq(2,length(rmse[[h]]),6)
    dizi_z = seq(3,length(rmse[[h]]),6)
    dizi = c(dizi_x, dizi_y, dizi_z)
    dizi = sort(dizi)
    rmse1 [[h]] = as.matrix(rmse [[h]] [dizi])
    rmse2 [[h]] = rmse[[h]][dizi_x]
    rmse2 [[h]] = cbind(rmse2[[h]],rmse[[h]][dizi_y],rmse[[h]][dizi_z])
    colnames(rmse2[[h]]) = c("RMSe_X","RMSe_Y","RMSe_Z")
    rownames(rmse2[[h]]) = c(pid[[h]])
    rmse2 [[h]] = rmse2 [[h]] [order(row.names(rmse2 [[h]]), decreasing = FALSE), ]
    rnms_rms [[h]] = row.names(rmse2[[h]])
    j = 1
    for (i in seq(1, length(rmse1 [[h]]), 3)){
      rmse1 [[h]] [i,1] = rmse2 [[h]] [rnms_rms [[h]] [j], "RMSe_X"]
      j = j+1
    }
    j = 1
    for (i in seq(2, length(rmse1 [[h]]), 3)){
      rmse1 [[h]] [i,1] = rmse2 [[h]] [rnms_rms [[h]] [j], "RMSe_Y"]
      j = j+1
    }
    j = 1
    for (i in seq(3, length(rmse1 [[h]]), 3)){
      rmse1 [[h]] [i,1] = rmse2 [[h]] [rnms_rms [[h]] [j], "RMSe_Z"]
      j = j+1
    }
    #print(rmse2[[h]])
  }
  
  #Processing the *.SNX files
    # First, required variables are defined
  sindex = list()
  ind1 = list()
  typ1 = list()
  cod1 = list()
  sindex1 = list()
  qvalues = list()
  para1 = list()
  para2 = list()
  para3 = list()
  para4 = list()
  para5 = list()
   
  # Station and coordinate indices are taken from whole data
  for(h in 1:length(variant2)){
    
    start1[h] = grep("+SOLUTION/ESTIMATE",variant2[[h]])[1] + 1
    finish1[h] = grep("-SOLUTION/ESTIMATE",variant2[[h]])[1] - 1
    sindex [[h]] = variant2[[h]][start1[h]:finish1[h]]
    
    start1[h] = regexpr("*INDEX",sindex[[h]][1])
    finish1[h] = start1[h] + nchar("*INDEX") -1
    
    ind1 [[h]] = matrix(as.numeric(substr(sindex[[h]][2],start1[h],finish1[h])))
    for (k in 3:length(sindex[[h]])){
      ind1 [[h]] = rbind(ind1 [[h]], matrix(as.numeric(substr(sindex[[h]][k],start1[h],finish1[h]))))
    }
    
    start1[h] = regexpr("TYPE__",sindex[[h]][1])
    finish1[h] = start1[h] + nchar("TYPE__")
    
    typ1 [[h]] = matrix(substr(sindex[[h]][2],start1[h],finish1[h]))
    for (k in 3:length(sindex[[h]])){
      typ1 [[h]] = rbind(typ1[[h]], matrix(substr(sindex[[h]][k],start1[h],finish1[h])))
    }
    typ1 [[h]] = matrix(str_extract(typ1[[h]],"[A-Z]+"))
    
    start1[h] = regexpr("CODE",sindex[[h]][1])
    finish1[h] = start1[h] + nchar("CODE")
    
    cod1 [[h]] = matrix(substr(sindex[[h]][2],start1[h],finish1[h]))
    for (k in 3:length(sindex[[h]])){
      cod1 [[h]] = rbind(cod1[[h]], matrix(substr(sindex[[h]][k],start1[h],finish1[h])))
    }
    cod1 [[h]] = matrix(str_extract(cod1[[h]],"[[A-Z]||[:digit:]][[A-Z]||[:digit:]][[A-Z]||[:digit:]][[A-Z]||[:digit:]]"))
    
    sindex1 [[h]] = cbind(ind1[[h]],typ1[[h]],cod1[[h]])
    
    #print(sindex1[[h]])
    
  }
  
  # Taking the elements of estimated covariance matrix from whole data
  for (h in 1:length(variant2)){
    
    start1[h] = grep("+SOLUTION/MATRIX_ESTIMATE L COVA",variant2[[h]])[1] + 1
    finish1[h] = grep("-SOLUTION/MATRIX_ESTIMATE L COVA",variant2[[h]])[1] - 1
    qvalues [[h]] =  variant2[[h]][start1[h]:finish1[h]]
    
    # Separating the elements regarding the estimated covariance matrix into columns
    start1[h] = regexpr("PARA1",qvalues[[h]][1])
    finish1[h] = nchar("PARA1") + start1[h]
    para1[[h]] = matrix(str_trim(substr(qvalues[[h]][2],start1[h],finish1[h])))
    for (k in 3:length(qvalues[[h]])){
      para1[[h]] = rbind(para1[[h]], matrix(str_trim(substr(qvalues[[h]][k],start1[h],finish1[h]))))
    }
    
    start1[h] = regexpr("PARA2",qvalues[[h]][1])
    finish1[h] = nchar("PARA2") + start1[h]
    para2 [[h]] = matrix(str_trim(substr(qvalues[[h]][2],start1[h],finish1[h])))
    for (k in 3:length(qvalues[[h]])){
      para2[[h]] = rbind(para2[[h]], matrix(str_trim(substr(qvalues[[h]][k],start1[h],finish1[h]))))
    }
    
    start1[h] = regexpr("[_][_][_][_][P][A][R][A][2][+][0]",qvalues[[h]][1])
    finish1[h] = nchar("____PARA2+0__________") + start1[h]
    para3 [[h]] = matrix(as.numeric(str_trim(substr(qvalues[[h]][2],start1[h],finish1[h]))))
    for (k in 3:length(qvalues[[h]])){
      para3[[h]] = rbind(para3[[h]], matrix(as.numeric(str_trim(substr(qvalues[[h]][k],start1[h],finish1[h])))))
    }
    
    start1[h] = regexpr("[_][_][_][_][P][A][R][A][2][+][1]",qvalues[[h]][1])
    finish1[h] = nchar("____PARA2+1__________") + start1[h]
    para4 [[h]] = matrix(as.numeric(str_trim(substr(qvalues[[h]][2],start1[h],finish1[h]))))
    for (k in 3:length(qvalues[[h]])){
      para4 [[h]] = rbind(para4[[h]], matrix(as.numeric(str_trim(substr(qvalues[[h]][k],start1[h],finish1[h])))))
    }
    
    start1[h] = regexpr("[_][_][_][_][P][A][R][A][2][+][2]",qvalues[[h]][1])
    finish1[h] = nchar("____PARA2+2__________") + start1[h]
    para5 [[h]] = matrix(as.numeric(str_trim(substr(qvalues[[h]][2],start1[h],finish1[h]))))
    for (k in 3:length(qvalues[[h]])){
      para5 [[h]] = rbind(para5[[h]], matrix(as.numeric(str_trim(substr(qvalues[[h]][k],start1[h],finish1[h])))))
    }
    
  }
  
  # Forming the estimated covariance matrix
  qmat = list()
  for (h in 1:length(variant2)){
    qmat[[h]] = matrix(para3[[h]][1,1],nrow = as.numeric(para1[[h]][(length(para1[[h]]) - 0):length(para1[[h]])]), ncol = as.numeric(para1[[h]][(length(para1[[h]]) - 0):length(para1[[h]])]))
    dimnames(qmat[[h]]) = list(1:para1[[h]][(length(para1[[h]]) - 0):length(para1[[h]])], 1:para1[[h]][(length(para1[[h]]) - 0):length(para1[[h]])])
    
    for(k in 1:length(para3[[h]])){
      qmat[[h]][para1[[h]][k,],para2[[h]][k,]] = para3[[h]][k,]
    }
    
    for(k in 1:length(para3[[h]])){
      qmat[[h]][para1[[h]][k,],matrix(as.numeric(para2[[h]])+1)[k,]] = para4[[h]][k,]
    }
    
    for(k in 1:length(para3[[h]])){
      qmat[[h]][para1[[h]][k,],matrix(as.numeric(para2[[h]])+2)[k,]] = para5[[h]][k,]
    }
    
    for(k in 2:length(para3[[h]])){
      qmat[[h]][para2[[h]][k,],para1[[h]][k,]] = para3[[h]][k,]
    }
    
    for(k in 1:length(para3[[h]])){
      if(is.na(para4[[h]][k,]) == FALSE){
        qmat[[h]][matrix(as.numeric(para2[[h]])+1)[k,], para1[[h]][k,]] = para4[[h]][k,]
      }
    }
    
    for(k in 1:length(para3[[h]])){
      if(is.na(para5[[h]][k,]) == FALSE){
        qmat[[h]][matrix(as.numeric(para2[[h]])+2)[k,],para1[[h]][k,]] = para5[[h]][k,]
      }
    }
  }
  
  #Alphabetical ordering of the elements of the estimated covariance matrix by station names
  rws = list()
  cls = list()
  qmat1 = list()
  rws1 = list()
  cls1 = list()
  for (h in 1:length(variant2)){
    rws[[h]] = paste(sindex1[[h]][,3], sindex1[[h]][,2], sep = "")
    cls[[h]] = paste(sindex1[[h]][,3], sindex1[[h]][,2], sep = "")
    colnames(qmat[[h]]) = cls[[h]]
    rownames(qmat[[h]]) = rws[[h]]
    qmat1[[h]] = matrix(rep(0,length(rws[[h]])*length(cls[[h]])), length(rws[[h]]), length(rws[[h]]))
    rws1[[h]] = sort(rws[[h]])
    cls1[[h]] = sort(cls[[h]])
    colnames(qmat1[[h]]) = cls1[[h]]
    rownames(qmat1[[h]]) = rws1[[h]]
    for (k in 1:length(rws1[[h]])){
      for (i in 1:length(cls1[[h]])){
        qmat1[[h]][rws1[[h]][k],cls1[[h]][i]] = qmat[[h]][rws1[[h]][k],cls1[[h]][i]]
      }
    }
  }
  
  #Exporting the dates of epochs
  dates = list()
  for(h in 1:length(variant)){
    dates [[h]] = variant[[h]][grep("Reference epoch:",variant[[h]])]
    dates [[h]] = readr::parse_date(str_extract(dates[h],"\\d\\d\\d\\d\\-\\d\\d\\-\\d\\d"), format = "%Y-%m-%d")
  }
  dates1 = apply( cbind( t(dates) ) , 1, unlist)
  orddates = cbind(sort(dates1, decreasing = FALSE))
  matcdates = cbind(match(orddates, dates1))
  
  #Matching the points between treated epochs
  both = list()
  onlyfirst = list()
  onlysecond = list()
  for (h in 2:length(variant)){
    for (t in 1:(length(variant)-1)){
      both [[t]] = pid[[which.min(dates)]][pid[[which.min(dates)]] %in% pid[[matcdates[h]]]] # in both, same as call: intersect(first, second)
      onlyfirst[[t]] = pid[[which.min(dates)]][!pid[[which.min(dates)]] %in% pid[[matcdates[h]]]] # only in 'first', same as: setdiff(first, second)
      onlysecond[[t]] = pid[[matcdates[h]]][!pid[[matcdates[h]]] %in% pid[[which.min(dates)]]] # only in 'second', same as: setdiff(second, first)
    }
  }
  
  #Carrying out IWST
  if (onlyfirst == "character(0)" && onlysecond == "character(0)"){
    
    #Variance ratio test
     #Calculating the pooled variance factors
     mortak = list()
     if ("TRUE" %in% duplicated(matcdates)){
       t = 1
       for (h in 1:length(variant)){
         if (h != which.min(dates)){
           mortak[[t]] = sqrt((Dof[which.min(dates)]*(rms[which.min(dates)])^2+Dof[h]*(rms[h])^2)/(Dof[which.min(dates)]+Dof[h]))
           t = t+1
         }
       }
     }else {
       t = 1
       for (h in 2:length(variant)){
         mortak[[t]] = sqrt((Dof[which.min(dates)]*(rms[which.min(dates)])^2+Dof[matcdates[h]]*(rms[matcdates[h]])^2)/(Dof[which.min(dates)]+Dof[matcdates[h]]))
         t = t+1
       }
     }
     
     #Calculating the pooled degree of freedoms
     dofort = list()
     if ("TRUE" %in% duplicated(matcdates)){
       t = 1
       for (h in 1:length(variant)){
         if (h != which.min(dates)){
           dofort[[t]] = (Dof[which.min(dates)]+Dof[h])
           t = t+1
         }
       }
     }else {
       t = 1
       for (h in 2:length(variant)){
         dofort[[t]] = (Dof[which.min(dates)]+Dof[matcdates[h]])
         t = t+1
       }
     }
     
     #Calculating the test statistics
     tval = list()
     if ("TRUE" %in% duplicated(matcdates)){
       t = 1
       for (h in 1:length(variant)){
         if (rms[which.min(dates)] > rms[h]){
           if (h != which.min(dates)){
             tval[[t]] = (rms[which.min(dates)]^2)/(rms[h]^2)
             t = t + 1
           }
         }else {
           if (h != which.min(dates)){
             tval[[t]] = (rms[h]^2)/(rms[which.min(dates)]^2)
             t = t + 1
           }
         }
       }
     }else {
       t = 1
       for (h in 2:length(variant)){
         if (rms[matcdates[1]] > rms[matcdates[h]]){
           tval[[t]] = (rms[matcdates[1]]^2)/(rms[matcdates[h]]^2)
           t = t + 1
         }else {
           tval[[t]] = (rms[matcdates[h]]^2)/(rms[matcdates[1]]^2)
           t = t + 1
         }
       }
     }
     
     tval1 = list()
     if ("TRUE" %in% duplicated(matcdates)){
       t = 1
       for (h in 1:length(rms)){
         if (abs(rms[which.min(dates)] - rms[h]) < 0.00051){
           if (h != which.min(dates)){
             tval1[[t]] = c("Variance ratio test passed")
             t = t + 1
           }
         }else {
           if (h != which.min(dates)){
             tval1[[t]] = c("Variance ratio test failed! Please carry out the adjustment computation again, after checking the adjustment parameters and stochastic model.")
             t = t + 1
           }
         }
       }
     }else {
       t = 1
       for (h in 2:length(rms)){
         if (abs(rms[matcdates[1]] - rms[matcdates[h]]) < 0.00051){
           tval1[[t]] = c("Variance ratio test passed")
           t = t + 1
         }else {
           tval1[[t]] = c("Variance ratio test failed! Please carry out the adjustment computation again, after checking the adjustment parameters and stochastic model.")
           t = t + 1
         }
       }
     }
     
     #Comparing test statistics to critical values
     cval = list()
     if ("TRUE" %in% duplicated(matcdates)){
       t = 1
       for (h in 1:length(variant)){
         if (h != which.min(dates)){
           cval[[t]] = qf(0.975,Dof[which.min(dates)],Dof[h])
           t = t + 1
         }
       }
     }else {
       t = 1
       for (h in 2:length(variant)){
         cval[[t]] = qf(0.975,Dof[matcdates[1]],Dof[matcdates[h]])
         t = t + 1
       }
     }
     
     pass = list()
     for (h in 1:length(cval)){
       if (as.numeric(cval[h]) > as.numeric(tval[h])){
         pass [[h]] = rbind(pass[h][!sapply(pass[h],is.null)], c("Variance ratio test passed"))
         cat(sprintf("Variance ratio test passed\n"))
         cat(sprintf("\n"))
       }else if(as.numeric(cval[h]) < as.numeric(tval[h])){
         pass [[h]] = rbind(pass[h][!sapply(pass[h],is.null)], tval1[[h]])
         cat(sprintf("%s\n",tval1[[h]]))
         cat(sprintf("\n"))
       }else {
         pass [[h]] = rbind(pass[h][!sapply(pass[h],is.null)], c("Variance ratio test failed! Please carry out the adjustment computation again, after checking the adjustment parameters and stochastic model."))
         cat(sprintf("Variance ratio test failed! Please carry out the adjustment computation again, after checking the adjustment parameters and stochastic model.\n"))
         cat(sprintf("\n"))
       }
     }
     
    #Calculating the displacement vectors
    d = list()
    if ("TRUE" %in% duplicated(matcdates)){
      t = 1
      for(h in 1:(length(v1))){
        if (h != which.min(dates)){
          d [[t]] = (v1[[h]] - v1[[which.min(dates)]])*1000
          t = t + 1
        }
      }
    }else {
      t = 1
      for(h in 2:(length(v1))){
        d [[t]] = (v1[[matcdates[h]]] - v1[[matcdates[1]]])*1000
        t = t + 1
      }
    }
    
    #Calculating the cofactor matrices of displacement vectors
    qd = list()
    if ("TRUE" %in% duplicated(matcdates)){
      t = 1
      for(h in 1:(length(qmat1))){
        if (h != which.min(dates)){
          qd [[t]] = (qmat1[[h]]/(rms[h]^2) + qmat1[[which.min(dates)]]/(rms[which.min(dates)]^2))
          t = t + 1
        }
      }
    }else {
      t = 1
      for(h in 2:(length(qmat1))){
        qd [[t]] = (qmat1[[matcdates[h]]]/(rms[matcdates[h]]^2) + qmat1[[matcdates[1]]]/(rms[matcdates[1]]^2))
        t = t + 1
      }
    }
    
    #Forming Gt matrices
    Gt = list()
    t = 1
    for(h in 1:(length(both))){
      Gt[[t]] = matrix(diag(3), 3, 3*length(both[[h]]))
      t = t + 1
    }
    
    #Forming identity matrices
    im = list()
    t = 1
    for(h in 1:(length(both))){
      im[[t]] = diag(3*length(both[[h]]))
      t = t + 1
    }
    
    #Calculating transformed displacements
    fark = list()
    W = list()
    S = list()
    dy = list()
    qy = list()
    #bilgi = list()
    
    t = 1
    for(h in 1:length(both)){
      W[[t]] = diag(3*length(both[[h]]))
      t = t + 1
      
      S [[h]] = im[[h]] - ((t(Gt[[h]])%*%(Ginv(Gt[[h]]%*%W[[h]]%*%t(Gt[[h]]))))%*%Gt[[h]]%*%W[[h]])
      d [[h]] = S[[h]] %*% d[[h]]
      qd [[h]] = S[[h]] %*% qd[[h]] %*% t(S[[h]])
      
    }
    
    for(h in 1:length(both)){
      for (k in 1:(3*length(both[[h]]))){
        if (abs(d[[h]][k]) < 0.01){
          W[[h]][k,k] = 0
        }else {
          W[[h]][k,k] = (1/(abs(d[[h]][k])))
        }
      }
    }
    
    ctest = list()
    status = list()
    it = matrix(rep(1),length(both))
    for(h in 1:length(both)){
      S [[h]] = im[[h]] - ((t(Gt[[h]])%*%(Ginv(Gt[[h]]%*%W[[h]]%*%t(Gt[[h]]))))%*%Gt[[h]]%*%W[[h]])
      dy [[h]] = S[[h]] %*% d[[h]]
      qy [[h]] = S[[h]] %*% qd[[h]] %*% t(S[[h]])
      fark[[h]] = abs(dy[[h]] - d[[h]])
      
      while ("FALSE" %in% (fark[[h]] < 0.1)){
        d[[h]] = dy[[h]]
        qd[[h]] = qy[[h]]
        for (k in 1:(3*length(both[[h]]))){
          if (abs(d[[h]][k]) < 0.01){
            W[[h]][k,k] = 0
          }else {
            W[[h]][k,k] = (1/(abs(d[[h]][k])))
          }
        }
        it[h,1] = it[h,1]+1
        S [[h]] = im[[h]] - ((t(Gt[[h]])%*%(Ginv(Gt[[h]]%*%W[[h]]%*%t(Gt[[h]]))))%*%Gt[[h]]%*%W[[h]])
        dy [[h]] = S[[h]] %*% d[[h]]
        qy [[h]] = S[[h]] %*% qd[[h]] %*% t(S[[h]])
        fark[[h]] = abs(dy[[h]] - d[[h]])
        if (it[h,1]>20) break
      }
      
      #Calculating the test statistics of coordinate components
      ctest[[h]] = matrix(((abs(dy[[h]][1,1])/1000)^2)/(qy[[h]][1,1]*(mortak[[h]]^2)),3*length(both[[h]]),1)
      for(t in 2:length(dy[[h]])){
        ctest[[h]][t,1] = ((abs(dy[[h]][t,1])/1000)^2)/(qy[[h]][t,t]*(mortak[[h]]^2))
      }
      
      #Determining the status of the points
      diz = seq(1,length(dy[[h]]),3)
      status[[h]] = na.exclude(matrix())
      for (t in 1:(length(dy[[h]])/3)){
        if("TRUE" %in% (qf(0.95,1,dofort[[h]]) < ctest[[h]][diz[t]]) | "TRUE" %in% (qf(0.95,1,dofort[[h]]) < ctest[[h]][diz[t]+1]) | "TRUE" %in% (qf(0.95,1,dofort[[h]]) < ctest[[h]][diz[t]+2])){
          status[[h]] = rbind(status[[h]],"unstable")
        }else{
          status[[h]] = rbind(status[[h]],"stable")
        }
      
      }
        
    }
    
    #Carrying out Final S transformation
    dfin = list()
    qdfin = list()
    Wfin = list()
    Sfin = list()
    
    t = 1
    for(h in 1:length(both)){
      Wfin[[t]] = diag(3*length(both[[h]]))
      t = t + 1
    }
    
    #Forming final weight matrix
    for(h in 1:length(both)){
      diz = seq(1,length(dy[[h]]),3)
      for(t in 1:length(status[[h]])){
        if (status[[h]][t] == "unstable"){
          for (k in diz[t]:(diz[t]+2)){
            Wfin[[h]][k,k] = 0
          }
        }else {
          for (k in diz[t]:(diz[t]+2)){
            Wfin[[h]][k,k] = 1
          }
        }
      }
    }
    
    #Calculating final displacements and cofactors
    for(h in 1:length(both)){
      Sfin [[h]] = im[[h]] - ((t(Gt[[h]])%*%(Ginv(Gt[[h]]%*%Wfin[[h]]%*%t(Gt[[h]]))))%*%Gt[[h]]%*%Wfin[[h]])
      dfin [[h]] = Sfin[[h]] %*% dy[[h]]
      qdfin [[h]] = Sfin[[h]] %*% qy[[h]] %*% t(Sfin[[h]])
    }
    
    #Splitting of qdfin into sub-matrices and dfin into sub-vectors
    krn = list()
    sbmt = list()
    krn1 = list()
    sbvt1 = list()
    for(h in 1:length(both)){
      krn[[h]] = kronecker(matrix(1:(dim(status[[h]])[1]^2), dim(status[[h]])[1], byrow = TRUE), matrix(1, 3, 3))
      sbmt[[h]] = lapply(split(qdfin[[h]], krn[[h]]), matrix, nr = 3)
      krn1[[h]] = kronecker(matrix(1:(dim(status[[h]])[1]), dim(status[[h]])[1], byrow = TRUE), matrix(1, 3, 1))
      sbvt1[[h]] = lapply(split(dfin[[h]], krn1[[h]]), matrix, nr = 3)
    }
    
    #Getting the 3 x 3 diagonal elements of qdfin
    diaqdfin = list()
    sbvt = list()
    for(h in 1:length(both)){
      dgnmbr = seq(1,(dim(status[[h]])[1]^2),(dim(status[[h]])[1]+1))
      diaqdfin[[h]] = na.exclude(matrix())
      for(k in 1:length(dgnmbr)){
        diaqdfin[[h]][k] = sbmt[[h]][dgnmbr[k]]
      }
      sbvt[[h]] = na.exclude(matrix())
      for(k in 1:dim(status[[h]])[1]){
        sbvt[[h]][k] = sbvt1[[h]][k]
      }
    }
    
    #calculating the absolute displacements
    dplc = list()
    dndedu = list()
    dplc1 = list()
    tist = list()
    statusfin = list()
    restab = list()
    dndedu1 = list()
    restab1 = list()
    for(h in 1:length(both)){
      dplc[[h]] = na.exclude(matrix())
      for(k in 1:dim(status[[h]])[1]){
        dplc[[h]] = rbind(dplc[[h]], sqrt(t(sbvt[[h]][[k]]) %*% sbvt[[h]][[k]]))
      }
      #calculating the single point test statistics
      tist[[h]] = na.exclude(matrix())
      for(k in 1:dim(status[[h]])[1]){
        tist[[h]] = rbind(tist[[h]], (t(sbvt[[h]][[k]]/1000) %*% Ginv(diaqdfin[[h]][[k]]) %*% (sbvt[[h]][[k]]/1000))/(3*mortak[[h]]^2))
      }
      #Forming the status vector of single point test
      statusfin[[h]] = na.exclude(matrix())
      for(k in 1:dim(status[[h]])[1]){
        if("TRUE" %in% (qf(0.95,3,dofort[[h]]) < tist[[h]][k])){
          statusfin[[h]] = rbind(statusfin[[h]],"unstable")
        }else{
          statusfin[[h]] = rbind(statusfin[[h]],"stable")
        }
      }
      #Calculating local geodetic displacements
      radius = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      aplat = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      aproxn = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      aproxh = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      newlat = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      newn = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      newh = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      newlat1 = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      lon = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      dif1 = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      it2 = matrix(rep(1),length(epochs[[which.min(dates)]][,1]),1)
      onon = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      ontw = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      onth = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      twon = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      twtw = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      twth = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      thon = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      thtw = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      thth = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      Rdon = matrix(rep(0),length(epochs[[which.min(dates)]][,1]),1)
      #GRS-80 ellipsoid parameters a = 6378137 m, b = 6356752.3141 m
      smaja = 6378137
      smina = 6356752.3141
      excen = ((smaja^2-smina^2)/smaja^2)
      for(k in 1:length(epochs[[which.min(dates)]][,1])){
        radius [k,1] = sqrt(epochs[[which.min(dates)]][k,1]^2+epochs[[which.min(dates)]][k,2]^2)
        aplat [k,1] = atan((epochs[[which.min(dates)]][k,3]/radius[k,1])*(1-excen)^-1)
        aproxn [k,1] = smaja^2/sqrt(smaja^2*cos(aplat[k,1])^2+smina^2*sin(aplat[k,1])^2)
        aproxh [k,1] = (radius[k,1]/cos(aplat[k,1]))-aproxn[k,1]
        newlat [k,1] = atan((epochs[[which.min(dates)]][k,3]/radius[k,1])*(1-((excen*aproxn[k,1])/(aproxn[k,1]+aproxh[k,1])))^-1)
        if (newlat[k,1] != aplat[k,1]){
          newn [k,1] = smaja^2/sqrt(smaja^2*cos(newlat[k,1])^2+smina^2*sin(newlat[k,1])^2)
          newh [k,1] = (radius[k,1]/cos(newlat[k,1]))-newn[k,1]
          newlat1 [k,1] = atan((epochs[[which.min(dates)]][k,3]/radius[k,1])*(1-((excen*newn[k,1])/(newn[k,1]+newh[k,1])))^-1)
          dif1 [k,1] = abs(newlat1 [k,1] - newlat [k,1])
          while (dif1 [k,1] > 0.000000000000001){
            newlat[k,1] = newlat1[k,1]
            newn [k,1] = smaja^2/sqrt(smaja^2*cos(newlat[k,1])^2+smina^2*sin(newlat[k,1])^2)
            newh [k,1] = (radius[k,1]/cos(newlat[k,1]))-newn[k,1]
            newlat1 [k,1] = atan((epochs[[which.min(dates)]][k,3]/radius[k,1])*(1-((excen*newn[k,1])/(newn[k,1]+newh[k,1])))^-1)
            dif1 [k,1] = abs(newlat1 [k,1] - newlat [k,1])
            it2[k,1] = it2[k,1]+1
          }
        }else{
          newlat1 [k,1] = aplat [k,1]
        }
        lon [k,1] = atan(epochs[[which.min(dates)]][k,2]/epochs[[which.min(dates)]][k,1])
        onon [k,1] = -sin(newlat1[k,1])*cos(lon[k,1])
        ontw [k,1] = -sin(lon[k,1])
        onth [k,1] = cos(newlat1[k,1])*cos(lon[k,1])
        twon [k,1] = -sin(newlat1[k,1])*sin(lon[k,1])
        twtw [k,1] = cos(lon[k,1])
        twth [k,1] = cos(newlat1[k,1])*sin(lon[k,1])
        thon [k,1] = cos(newlat1[k,1])
        thtw [k,1] = 0
        thth [k,1] = sin(newlat1[k,1])
        for (i in 1:length(epochs[[which.min(dates)]][,1])){
          Rdon [i] = list(matrix(rep(0),3,3))
        }
        for (j in 1: length(epochs[[which.min(dates)]][,1])){
          Rdon[[j]][1,1] = onon [j,1]
          Rdon[[j]][1,2] = ontw [j,1]
          Rdon[[j]][1,3] = onth [j,1]
          Rdon[[j]][2,1] = twon [j,1]
          Rdon[[j]][2,2] = twtw [j,1]
          Rdon[[j]][2,3] = twth [j,1]
          Rdon[[j]][3,1] = thon [j,1]
          Rdon[[j]][3,2] = thtw [j,1]
          Rdon[[j]][3,3] = thth [j,1]
        }
      }
      geocoor1 = list()
      for(k in 1:length(epochs[[which.min(dates)]][,1])){
        geocoor1[[k]] =  t(Rdon[[k]]) %*% sbvt1[[h]][[k]]
      }
      dndedu[[h]] = geocoor1
      dndedu1[[h]] = apply( cbind( t(dndedu[[h]]) ) , 1, unlist)
      dplc1[[h]] = na.exclude(matrix())
      for(k in 1:length(epochs[[which.min(dates)]][,1])){
        dplc1[[h]] = rbind(dplc1[[h]], sqrt(t(dndedu[[h]][[k]]) %*% dndedu[[h]][[k]]))
      }
      
      #Forming the results tables
      restab[[h]] = matrix(both[[h]],length(status[[h]]),1)
      restab[[h]] = cbind(restab[[h]],formatC(matrix(dfin[[h]][diz],length(status[[h]]),1),digits = 4,format = "f"))
      restab[[h]] = cbind(restab[[h]],formatC(matrix(dfin[[h]][diz+1],length(status[[h]]),1),digits = 4,format = "f"))
      restab[[h]] = cbind(restab[[h]],formatC(matrix(dfin[[h]][diz+2],length(status[[h]]),1),digits = 4,format = "f"))
      restab[[h]] = cbind(restab[[h]],formatC(dplc[[h]],digits = 4,format = "f"))
      restab[[h]] = cbind(restab[[h]],formatC(tist[[h]],digits = 4,format = "f"))
      restab[[h]] = cbind(restab[[h]],statusfin[[h]])
      rownames(restab[[h]]) = c(both[[h]])
      colnames(restab[[h]]) = c("Station", "DX (mm)", "DY (mm)","DZ (mm)","Disp Vector (mm)","Test Value","Status")
      print(restab[[h]], digits = 4, na.print = "",zero.print = "0")
      cat(sprintf("\n"))
      
      restab1[[h]] = matrix(both[[h]],length(status[[h]]),1)
      restab1[[h]] = cbind(restab1[[h]],formatC(matrix(dndedu1[[h]][diz],length(status[[h]]),1),digits = 4,format = "f"))
      restab1[[h]] = cbind(restab1[[h]],formatC(matrix(dndedu1[[h]][diz+1],length(status[[h]]),1),digits = 4,format = "f"))
      restab1[[h]] = cbind(restab1[[h]],formatC(matrix(dndedu1[[h]][diz+2],length(status[[h]]),1),digits = 4,format = "f"))
      restab1[[h]] = cbind(restab1[[h]],formatC(dplc1[[h]],digits = 4,format = "f"))
      rownames(restab1[[h]]) = c(both[[h]])
      colnames(restab1[[h]]) = c("Station", "Dn (mm)", "De (mm)","Du (mm)","Disp Vector (mm)")
      print(restab1[[h]], digits = 4, na.print = "",zero.print = "0")
      cat(sprintf("\n"))
      
      #Forming the Dof lists
      dflist = na.exclude(matrix())
      if ("TRUE" %in% duplicated(matcdates)){
        for (p in 1:length(epochs)){
          dflist = rbind(dflist, c(sprintf("Degree of Freedom of Epoch on %s : %i", format.Date (dates[[p]],format = "%d-%m-%Y"), Dof[p])))
        }
      }else {
        for (p in 1:length(epochs)){
          dflist = rbind(dflist, c(sprintf("Degree of Freedom of Epoch on %s : %i", format.Date (dates[[matcdates[p]]],format = "%d-%m-%Y"), Dof[matcdates[p]])))
        }
      }
      
      #Forming the rms lists
      rmslist = na.exclude(matrix())
      if ("TRUE" %in% duplicated(matcdates)){
        for (p in 1:length(rms)){
          rmslist = rbind(rmslist, c(sprintf("A Posteriori RMS of Unit Weight of Epoch on %s : %s m", format.Date (dates[[p]],format = "%d-%m-%Y"), format(rms[p],digits = 5))))
        }
      }else {
        for (p in 1:length(rms)){
          rmslist = rbind(rmslist, c(sprintf("A Posteriori RMS of Unit Weight of Epoch on %s : %s m", format.Date (dates[[matcdates[p]]],format = "%d-%m-%Y"), format(rms[matcdates[p]],digits = 5))))
        }
      }
    }
    
     #Forming the pooled variance factor lists
     pvlist = na.exclude(matrix())
     tresult = na.exclude(matrix())
     if ("TRUE" %in% duplicated(matcdates)){
       for (p in 1:length(mortak)){
         pvlist = rbind(pvlist, c(sprintf("Pooled Variance Factor on %s : %s", format.Date (dates[[p+1]],format = "%d-%m-%Y"), format(round((as.numeric(mortak[[p]]^2)*1000000), 2), nsmall=2))))
         tresult = rbind(tresult,c(paste(pass[[p]],c(sprintf("for %s",format.Date (dates[[p+1]],format = "%d-%m-%Y"))))))
       }
     }else {
       for (p in 1:length(mortak)){
         pvlist = rbind(pvlist, c(sprintf("Pooled Variance Factor on %s : %s", format.Date (dates[[matcdates[p+1]]], format = "%d-%m-%Y"), format(round((as.numeric(mortak[[p]]^2)*1000000), 2), nsmall=2))))
         tresult = rbind(tresult,c(paste(pass[[p]],c(sprintf("for %s",format.Date (dates[[matcdates[p+1]]], format = "%d-%m-%Y"))))))
     }
    
     Sys.setenv(R_GSCMD="C:/Program Files/gs/gs9.53.3/bin/gswin64.exe")
     rmarkdown::render("Reporting.Rmd","pdf_document")
     rmarkdown::render("Reporting.Rmd","html_document")
     rstudioapi::showDialog(title = "Information",
                            message = "The analyses were successfully <b>completed.</b>")
     }
     
  }else {
    cat(sprintf("Points could not be matched!\n %s points are missing in the first period and %s points in the second period.\n Please carry out the adjustment computation again and note that all files have the same points.", onlysecond, onlyfirst))
    rstudioapi::showDialog(title = "Warning!",
                           message = sprintf("Points could not be <b>matched!</b>\n <b>%s</b> points are missing in the first period and <b>%s</b> points in the second period.\n Please carry out the adjustment computation again and note that all files have the same points.", onlysecond, onlyfirst))
  }
  
}else {
  cat(sprintf("The number of selected .OUT and .SNX files have to be equal !!! "))
  rstudioapi::showDialog(title = "Warning!",
                         message = "The number of selected .OUT and .SNX files have to be <b>equal!</b>")
}
