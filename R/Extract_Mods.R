Extract_Mods = function(network_file, sem_score, neighbours_file, clusteringcoeff, CCT = 0.5, SST=0.5, min_clus_genes = 3) {
  network_Cnt =  network_file #network file
  semsc_Cnt = sem_score
  neighbours_Cnt = neighbours_file
  cc_Cnt = clusteringcoeff
  cc_Cnt = filter(cc_Cnt, ClusteringCoeff > CCT) #took care of threshold part for clustering coeffecient
  min_clus_genes = 3

  Modules = data.frame("ModuleNo" = numeric(), "Genes" = character())
  module_count = 1

  clus_coef = merge(cc_Cnt,neighbours_Cnt,by = "gene_id")
  clus_coef = arrange(clus_coef, desc(ClusteringCoeff))

  #module creation starts
  q =1


  while (nrow(clus_coef) > 4) {

    seed_gene = as.character(clus_coef$gene_id[1])
    temp_cluster = c(seed_gene)
    seed_neighbours = unlist(strsplit(as.character(clus_coef$Neighbours[1]),", "))
    counted_neighbours = c(seed_gene)
    prev_length = 1
    print("came here too")

    while (length(counted_neighbours) <= length(seed_neighbours)) {

      k = 1

      if (nrow(clus_coef) < 4) {
        break
      }

      repeat {

        i = seed_neighbours[k]


        if (is.na(i) || (nrow(clus_coef) < 4)) {
          break
        }

        if (!(i %in% counted_neighbours)) {

          for (x in temp_cluster) {

            if (class(try(sem_score[trimws(x),trimws(i)],silent = TRUE)) == "numeric") {

              print(sem_score[trimws(x),trimws(i)])

              if ((x != i) && (sem_score[trimws(x),trimws(i)] > SST)) {

                temp_cluster = c(temp_cluster, i)
                new_neighbours = unlist(strsplit(as.character(clus_coef$Neighbours[which(clus_coef$gene_id == as.character(i))]),","))
                seed_neighbours = c(seed_neighbours, new_neighbours)
                index = which(as.character(clus_coef$gene_id) == as.character(trimws(i)))
                clus_coef = clus_coef[-index, ]
                print("reduced x in i")
                break

              }

            }

          }

          counted_neighbours = c(counted_neighbours, i)

        }

        if (length(counted_neighbours) > length(seed_neighbours)) {

          clus_coef = clus_coef[-1, ]
          print("reduced vm")
          break

        }

        if (length(temp_cluster) > prev_length) {

          k = 1
          prev_length = length(temp_cluster)

        }
        else if (length(temp_cluster) == prev_length) {

          k = k + 1

        }

      }

    }

    if (length(temp_cluster) >= 3) {

      Modules = rbind(Modules,data.frame("ModuleNo" = module_count, "Genes" = toString(temp_cluster,sep = ",")))
      module_count = module_count + 1
      print("yahoo!!")
    }
    print(nrow(Modules))
  }



  return(Modules)

}
