################## The label correcting algorithm ####################
##### Copyright (C) 2016.11.16 Jaeyoung Jung all rights reserved. ####

rm(list=ls(all=TRUE))
gc(reset=TRUE)

### Convert the network information into a matrix - Chapter 5 pp.125 ~ 129
# Path
OD <- as.data.frame(matrix(c(1, 2, 2, 3, 3, 6, 1, 4, 4, 5,
                             5, 6, 1, 5, 5, 2, 2, 6, 2, 5), ncol = 2, byrow = TRUE))
colnames(OD) <- c("From", "To")
OD <- OD[order(OD$From, OD$To), ]
OD$Travel_time <- c(6, 3, 2, 2, 2, 1, 3, 1, 3, 5)
rownames(OD) <- c(1:10)

# Initial setting
Label_list <- rep(Inf, 6)
Predecessor_list <- rep(0, 6)
Label_list[1] <- 0
Sequence_list <- NULL
Sequence_list[1] <- 1
Origin <- 0
Destination <- 6

# Min path algorithm - Label correcting method
while (TRUE) {
  arr <- which(OD$From == Sequence_list[1])

  for (j in 1:length(arr)) {
    if (j < length(arr)) {
        if ((as.numeric(Label_list[OD[arr[j], 1]]) + OD[arr[j], 3]) < Label_list[OD[arr[j], 2]]) {
            Label_list[OD[arr[j], 2]] <- Label_list[OD[arr[j], 1]] + OD[arr[j], 3]
            Predecessor_list[OD[arr[j], 2]] <- OD[arr[j], 1]
            Sequence_list <- c(Sequence_list, OD[arr[j], 2])
        } else {
            Sequence_list <- Sequence_list[Sequence_list != 6] # Delete destination
        }
    } else {
        if ((as.numeric(Label_list[OD[arr[j], 1]]) + OD[arr[j], 3]) < Label_list[OD[arr[j], 2]]) {
            Label_list[OD[arr[j], 2]] <- Label_list[OD[arr[j], 1]] + OD[arr[j], 3]
            Predecessor_list[OD[arr[j], 2]] <- OD[arr[j], 1]
            Sequence_list <- c(Sequence_list, OD[arr[j], 2])
            Sequence_list <- Sequence_list[Sequence_list != OD[arr[j], 1]]
            Sequence_list <- Sequence_list[Sequence_list != 6] # Delete destination
        } else {
            Sequence_list <- Sequence_list[Sequence_list != OD[arr[j], 1]]
            Sequence_list <- Sequence_list[Sequence_list != 6] # Delete destination
        }
    }
  }
  print("EEEEEEEEEEEEEEEEEEEEEE")
  print(unlist(Label_list))
  print(unlist(Predecessor_list))
  print(unlist(Sequence_list))
  print("EEEEEEEEEEEEEEEEEEEEEE")
  Sys.sleep(3)
  if (length(Sequence_list) == 0) break
}

# Show the minimum travel time path
print(paste(Destination, " <- ", Predecessor_list[Destination], " <- ", Predecessor_list[Predecessor_list[Destination]],
            " <- ", Predecessor_list[Predecessor_list[Predecessor_list[Destination]]], " / ",
            "Minimum travel time is ", Label_list[Destination], sep = ""))

