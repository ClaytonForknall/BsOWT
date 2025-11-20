#' Generate B-spline wavelet transform matrix
#'
#' This function generates the B-spline wavelet transform matrix. It requires
#' the user to specify a quantitative explanatory variable, along with the type
#' of B-spline wavelet (one of `haar`, `linear` or `cubic`).
#'
#' @param x Vector corresponding to a quantitative explanatory variable
#' @param type Character string specifying the type of B-spline wavelet, being one of `haar`, `linear` or `cubic`
#' @return List of output. `W.tilde` contains the wavelet transform matrix
#' @export
BsOWT <- function(x = NULL, type = NULL){

  # type = "cubic"
  # x= seq(1, 10, 1)

  ##Check if the type is one of "haar", "linear", "cubic"
  if(!(type %in% c("haar","linear","cubic"))){
    stop("Option 'type' needs to be specified as 'haar','linear' or 'cubic'. No other types of B-spline are currently supported")
  }

  ##Check if x is numeric
  if(!(is.numeric(x))){
    stop("x has to be a numeric vector")
  }

  ##First set up a few variables which will be used throughout
  n <- length(x)
  J <- ceiling(log(n, base=2))
  orthog = TRUE

  ##define ptilde, which is the order of the basis function selected
  if(type == "haar"){
    ptilde <- 1
  } else if(type == "linear"){
    ptilde <- 2
  } else {
    ptilde <- 4
  }

  ##Now I'm going to put together a big loop, where in each iteration, all the elements of each level are calculated
  ##j here corresponds to a level of the wavelet transform and ranges from 0 to J (where 0 is coarsest and J is finest)

  ##Define things needed inside the loop and initialise
  ##J+1 index used to avoid the need for a 0 list index
  ##rather than specifying the below as simply list() without a length, I want to initialise each to have a length equal to J+1
  x_j <- vector(mode = "list", length = J+1)
  n_j <- vector(mode = "list", length = J+1)
  ne_j <- vector(mode = "list", length = J+1)
  no_j <- vector(mode = "list", length = J+1)
  P_j.list <- vector(mode = "list", length = J+1)
  D_j.list <- vector(mode = "list", length = J+1)
  U_j.list <- vector(mode = "list", length = J+1)
  H_j.list <- vector(mode = "list", length = J+1)
  G_j.list <- vector(mode = "list", length = J+1)
  H.tilde_j.list <- vector(mode = "list", length = J+1)
  H.tilde0_j.list <- vector(mode = "list", length = J+1)
  G.tilde_j.list <- vector(mode = "list", length = J+1)
  V_j.list <- vector(mode = "list", length = J+1)
  V.tilde_j.list <- vector(mode = "list", length = J+1)
  M_j.list <- vector(mode = "list", length = J+1)
  O_j.list <- vector(mode = "list", length = J+1)
  upsilon_j.list <- vector(mode = "list", length = J+1)
  sol.upsilon_j.list <- vector(mode = "list", length = J+1)
  W.tilde_j.disjoint <- vector(mode = "list", length = J+1)
  W.tilde_j <- vector(mode = "list", length = J+1)
  W.tilde_j.expand <- vector(mode = "list", length = J+1)
  temp.call <- vector(mode = "list", length = J+1)
  temp.expr <- NULL

  for(j in (J+1):2){

    ##initialise some variables, when j == J+1
    if(j == (J+1)){
      x_j[[j]] <- x
      n_j[[j]] <- n
    }

    ##first thing that I will do is build the even and odd vectors - these are used to split the data into even and odd observations
    ##much of what follows is built around these vectors
    ##although these labels seem back to front, remember x[1] is actually x_0
    even.index <- seq(1, n_j[[j]], by=2)
    odd.index <- seq(2, n_j[[j]], by=2)

    ##record the length of each vector
    ne_j[[j]] <- length(even.index)
    no_j[[j]] <- length(odd.index)

    ##in order to handle the expansion of the range of x's more easily, I'm going to set up a labelled vector
    xx_j <- c(rep(x_j[[j]][1], ptilde-1),  x_j[[j]], rep(x_j[[j]][n_j[[j]]], ptilde))
    names(xx_j) <- seq((even.index[1] - ptilde), (n_j[[j]] + ptilde -1), by = 1)

    ##now, let's form the different matrices needed in the wavelet transform
    ##some of these matrices will only have non-binary elements for the cubic B-spline
    ##the four matrices to be defined are:
    ##U - zero matrix for haar and linear, non-binary for cubic
    ##D - only non-binary for cubic
    ##P - non-binary for linear and cubic
    ##upsilon - non-binary for all types

    ################################
    ##U matrix
    ################################
    if(type == "cubic"){
      ##U matrix will have non-binary elements if type is cubic
      if(j==2){
        U_j.list[[j]] <- Matrix::Matrix(-1*(xx_j[names(xx_j) %in% (odd.index-1-1-1)] - xx_j[names(xx_j) %in% (odd.index-1-1-2)])/(xx_j[names(xx_j) %in% (odd.index-1-1+2)] - xx_j[names(xx_j) %in% (odd.index-1-1-2)]),
                                nrow=1, ncol=1)
      } else{
        diag.U.mat <- -1*(xx_j[names(xx_j) %in% (odd.index-1-1-1)] - xx_j[names(xx_j) %in% (odd.index-1-1-2)])/(xx_j[names(xx_j) %in% (odd.index-1-1+2)] - xx_j[names(xx_j) %in% (odd.index-1-1-2)])
        U_j.list[[j]] <- diag(diag.U.mat, nrow=ne_j[[j]], ncol=no_j[[j]])
      }
      ##for the cubic B-spline, we also have one off-diagonal element (lower triangular)
      off.diag.U.mat <- -1*(xx_j[names(xx_j) %in% (odd.index-1-1+2)] - xx_j[names(xx_j) %in% (odd.index-1-1+1)])/(xx_j[names(xx_j) %in% (odd.index-1-1+2)] - xx_j[names(xx_j) %in% (odd.index-1-1-2)])
      ##setting up a way to insert these elements - this will likely seem like overkill, but need to do this for when the matrix dimensions are different
      if(length(off.diag.U.mat) == 1){
        off.diag.U.mat <- 0
      } else if(ne_j[[j]] != no_j[[j]]){
        off.diag.U.mat <- c(off.diag.U.mat[2:length(off.diag.U.mat)],0)
      } else{
        off.diag.U.mat <- off.diag.U.mat[2:length(off.diag.U.mat)]
      }
      suppressWarnings(U_j.list[[j]][which(row(U_j.list[[j]]) == (col(U_j.list[[j]])+1))] <- off.diag.U.mat)
    } else{
      ##U matrix is just composed of zeroes if type is haar or linear
      U_j.list[[j]] <- Matrix::Matrix(0, nrow = ne_j[[j]], ncol = no_j[[j]], sparse = TRUE)
    }

    ################################
    ##D matrix
    ################################
    if(type == "cubic"){
      ##D matrix will have non-binary elements if type is cubic
      if(j==2){
        D_j.list[[j]] <- Matrix::Matrix((xx_j[names(xx_j) %in% (even.index-1+1)] - xx_j[names(xx_j) %in% (even.index-1-1)])/(xx_j[names(xx_j) %in% (even.index-1+2)] - xx_j[names(xx_j) %in% (even.index-1-2)]),
                                nrow=1, ncol=1)
      } else{
        D_j.list[[j]] <- Matrix::Diagonal(x = (xx_j[names(xx_j) %in% (even.index-1+1)] - xx_j[names(xx_j) %in% (even.index-1-1)])/(xx_j[names(xx_j) %in% (even.index-1+2)] - xx_j[names(xx_j) %in% (even.index-1-2)]))
      }
    } else{
      ##D matrix is just an identity if type is haar or linear
      D_j.list[[j]] <- Matrix::Diagonal(n = ne_j[[j]])
    }

    ################################
    ##P matrix
    ################################
    if(type == "cubic"){
      ##P matrix will have non-binary elements if type is cubic
      if(j==2){
        P_j.list[[j]] <- Matrix::Matrix((xx_j[names(xx_j) %in% (odd.index-1-1+4)] - xx_j[names(xx_j) %in% (odd.index-1-1+1)])/(xx_j[names(xx_j) %in% (odd.index-1-1+4)] - xx_j[names(xx_j) %in% (odd.index-1-1-2)]),
                                nrow=1, ncol=1)
      }else{
        diag.P.mat <- (xx_j[names(xx_j) %in% (odd.index-1-1+4)] - xx_j[names(xx_j) %in% (odd.index-1-1+1)])/(xx_j[names(xx_j) %in% (odd.index-1-1+4)] - xx_j[names(xx_j) %in% (odd.index-1-1-2)])
        P_j.list[[j]] <- Matrix::Matrix(diag(diag.P.mat, nrow=no_j[[j]], ncol=ne_j[[j]]), sparse = TRUE)
      }
      ##for the cubic B-spline, we also have one off-diagonal element (upper triangular)
      off.diag.P.mat <- (xx_j[names(xx_j) %in% (odd.index-1-1+1)] - xx_j[names(xx_j) %in% (odd.index-1-1-2)])/(xx_j[names(xx_j) %in% (odd.index-1-1+4)] - xx_j[names(xx_j) %in% (odd.index-1-1-2)])
      suppressWarnings(P_j.list[[j]][which(row(P_j.list[[j]]) == (col(P_j.list[[j]])-1))] <- off.diag.P.mat)

    } else if(type == "linear"){
      ##P matrix will have non-binary elements if type is linear
      if(j==2){
        P_j.list[[j]] <- Matrix::Matrix((xx_j[names(xx_j) %in% (odd.index-1-1+2)] - xx_j[names(xx_j) %in% (odd.index-1-1+1)])/(xx_j[names(xx_j) %in% (odd.index-1-1+2)] - xx_j[names(xx_j) %in% (odd.index-1-1)]),
                                nrow=1, ncol=1)
      }else{
        diag.P.mat <- (xx_j[names(xx_j) %in% (odd.index-1-1+2)] - xx_j[names(xx_j) %in% (odd.index-1-1+1)])/(xx_j[names(xx_j) %in% (odd.index-1-1+2)] - xx_j[names(xx_j) %in% (odd.index-1-1)])
        P_j.list[[j]] <- Matrix::Matrix(diag(diag.P.mat, nrow=no_j[[j]], ncol=ne_j[[j]]), sparse = TRUE)
      }
      ##for the linear B-spline, we also have one off-diagonal element (upper triangular)
      off.diag.P.mat <- (xx_j[names(xx_j) %in% (odd.index-1-1+1)] - xx_j[names(xx_j) %in% (odd.index-1-1)])/(xx_j[names(xx_j) %in% (odd.index-1-1+2)] - xx_j[names(xx_j) %in% (odd.index-1-1)])
      suppressWarnings(P_j.list[[j]][which(row(P_j.list[[j]]) == (col(P_j.list[[j]])-1))] <- off.diag.P.mat)

    } else{
      ##P matrix is just an identity if type is haar
      P_j.list[[j]] <- Matrix::Matrix(diag(nrow = no_j[[j]], ncol = ne_j[[j]]), sparse = TRUE)
    }

    ##before forming the upsilon matrix, there are some other matrices that need to be calculated first
    ##these are the HG_j matices

    ##next use equation 45 from Jansen 2016 (actually implemented using Eq 4.51 from Jansen 2022 textbook) to obtain the H and G matrices in augmented and disjoint form
    eq.45.mat <- (rbind(cbind(Matrix::Diagonal(n = ne_j[[j]]), -1*U_j.list[[j]]),
                        cbind(Matrix::Matrix(0, ncol = ne_j[[j]], nrow = no_j[[j]]), Matrix::Diagonal(n = no_j[[j]]))) %*%
                    rbind(cbind(D_j.list[[j]], Matrix::Matrix(0, nrow = ne_j[[j]], ncol = no_j[[j]])),
                          cbind(Matrix::Matrix(0, nrow = no_j[[j]], ncol = ne_j[[j]]), Matrix::Diagonal(n = no_j[[j]]))) %*%
                    rbind(cbind(Matrix::Diagonal(n = ne_j[[j]]), Matrix::Matrix(0, nrow = ne_j[[j]], ncol = no_j[[j]])),
                          cbind(P_j.list[[j]], Matrix::Diagonal(n = no_j[[j]]))))

    ##store a temporary version of H_j and G_j from this matrix
    temp.H_j <- suppressWarnings(Matrix::Matrix(eq.45.mat[, 1:ne_j[[j]]], ncol=ne_j[[j]]))
    temp.G_j <- suppressWarnings(Matrix::Matrix(eq.45.mat[, (ne_j[[j]]+1):n_j[[j]]], ncol=no_j[[j]]))
    ##the above are in disjoint form - need to convert them to sequential
    H_j.list[[j]] <- Matrix::Matrix(0, nrow=nrow(temp.H_j), ncol=ncol(temp.H_j))
    H_j.list[[j]][even.index,] <- temp.H_j[(1:ne_j[[j]]),]
    H_j.list[[j]][odd.index,] <- temp.H_j[(ne_j[[j]]+1):(n_j[[j]]),]
    ##now do G_j
    G_j.list[[j]] <- Matrix::Matrix(0, nrow=nrow(temp.G_j), ncol=ncol(temp.G_j))
    G_j.list[[j]][even.index,] <- temp.G_j[(1:ne_j[[j]]),]
    G_j.list[[j]][odd.index,] <- temp.G_j[(ne_j[[j]]+1):(n_j[[j]]),]

    ##also store H.tilde0_j and G.tilde_j - same deal as above. Will need to convert them from disjoint to sequential
    eq.45.mat.inv <- Matrix::solve(eq.45.mat)
    temp.H.tilde0_j <- suppressWarnings(Matrix::t(Matrix::Matrix(eq.45.mat.inv[1:ne_j[[j]],], nrow=ne_j[[j]])))
    temp.G.tilde_j <- suppressWarnings(Matrix::t(Matrix::Matrix(eq.45.mat.inv[(ne_j[[j]]+1):n_j[[j]],], nrow=no_j[[j]])))
    ##convert to sequential form
    H.tilde0_j.list[[j]] <- Matrix::Matrix(0, nrow=nrow(temp.H.tilde0_j), ncol=ncol(temp.H.tilde0_j))
    H.tilde0_j.list[[j]][even.index,] <- temp.H.tilde0_j[(1:ne_j[[j]]),]
    H.tilde0_j.list[[j]][odd.index,] <- temp.H.tilde0_j[(ne_j[[j]]+1):(n_j[[j]]),]
    G.tilde_j.list[[j]] <- Matrix::Matrix(0, nrow=nrow(temp.G.tilde_j), ncol=ncol(temp.G.tilde_j))
    G.tilde_j.list[[j]][even.index,] <- temp.G.tilde_j[(1:ne_j[[j]]),]
    G.tilde_j.list[[j]][odd.index,] <- temp.G.tilde_j[(ne_j[[j]]+1):(n_j[[j]]),]

    ##calculate V_j and V.tilde_j
    ##initialise V.tilde_j beyond the range of j to be the identity matrix
    ##note that V.tilde_j always operates at a scale j+1 finer than the current scale j

    if(j==J+1){
      #V_j.list[[j+1]] <- Matrix::Diagonal(n = n_j[[j]])
      V.tilde_j.list[[j+1]] <- Matrix::Diagonal(n = n_j[[j]])

      V_j.list[[j]] <- H_j.list[[j]]
      #V.tilde_j.list[[j]] <- t(H.tilde_j.list[[j]])

    } else{
      ##if j != J+1, then we need to form an expression that will be evaluated to calculate V_j and V.tilde_j
      ##deal with V_j first
      for(i in (J+1):j){
        if(i == J+1){
          temp.V_j.expr <- parse(text=paste0("H_j.list[[",i,"]]"))
        } else{
          temp.V_j.expr <- paste(temp.V_j.expr, parse(text=paste0("H_j.list[[",i,"]]")), sep="%*%")
        }
      }
      ##now sort V.tilde_j, remembering that this matrix operates at a scale behind V_j
      for(i in (J+1):(j+1)){
        if(i == J+1){
          temp.V.tilde_j.expr <- parse(text=paste0("Matrix::t(H.tilde_j.list[[",i,"]])"))
        } else{
          temp.V.tilde_j.expr <- paste(parse(text=paste0("Matrix::t(H.tilde_j.list[[",i,"]])")), temp.V.tilde_j.expr, sep="%*%")
        }
      }

      ##now form V_j and V.tilde_j
      V_j.list[[j]] <- eval(parse(text=temp.V_j.expr))
      V.tilde_j.list[[j+1]] <- eval(parse(text=temp.V.tilde_j.expr))
    }

    ##now calculate the moments of the spline function
    ##this calculation only has to be done once, when j==J+1
    ##however, the moments vary for the different spline functions
    if(j==J+1){
      if(type == "haar"){
        ##haar moments, obtained from Mathematica
        ##we're assuming there are one vanishing moment for the haar wavelet
        m0.vec <- (xx_j[names(xx_j) %in% (1:n_j[[j]]-1+1)] - xx_j[names(xx_j) %in% (1:n_j[[j]]-1)])
        M_j.list[[j+1]] <- matrix(m0.vec, nrow=length(m0.vec))

      } else if(type == "linear"){
        ##linear moments, obtained from Mathematica
        ##we're assuming there are two vanishing moments
        m0.vec <- (xx_j[names(xx_j) %in% (1:n_j[[j]]-1+1)] - xx_j[names(xx_j) %in% (1:n_j[[j]]-1-1)])/2
        m1.vec <- (xx_j[names(xx_j) %in% (1:n_j[[j]]-1+1)] - xx_j[names(xx_j) %in% (1:n_j[[j]]-1-1)])*(xx_j[names(xx_j) %in% (1:n_j[[j]]-1-1)] + xx_j[names(xx_j) %in% (1:n_j[[j]]-1)] +
                                                                                                         xx_j[names(xx_j) %in% (1:n_j[[j]]-1+1)])/6
        M_j.list[[j+1]] <- cbind(m0.vec, m1.vec)

      }else{
        ##cubic moments, obtained from Mathematica
        ##we're assuming there are two vanishing moments
        m0.vec <- (xx_j[names(xx_j) %in% (1:n_j[[j]]-1+2)] - xx_j[names(xx_j) %in% (1:n_j[[j]]-1-2)])/4
        m1.vec <- (xx_j[names(xx_j) %in% (1:n_j[[j]]-1+2)] - xx_j[names(xx_j) %in% (1:n_j[[j]]-1-2)])*(xx_j[names(xx_j) %in% (1:n_j[[j]]-1-2)] + xx_j[names(xx_j) %in% (1:n_j[[j]]-1-1)] +
                                                                                                         xx_j[names(xx_j) %in% (1:n_j[[j]]-1)] + xx_j[names(xx_j) %in% (1:n_j[[j]]-1+1)] +
                                                                                                         xx_j[names(xx_j) %in% (1:n_j[[j]]-1+2)])/20
        M_j.list[[j+1]] <- cbind(m0.vec, m1.vec)
      }
    }

    ##now define the M_j and O_j matrices at the current level j
    M_j.list[[j]] <- Matrix::Matrix(Matrix::t(H_j.list[[j]]) %*% M_j.list[[j+1]], sparse = TRUE)
    O_j.list[[j]] <- Matrix::Matrix(Matrix::t(G_j.list[[j]]) %*% M_j.list[[j+1]], sparse = TRUE)

    #####################################
    #####################################
    ### ORTHOGONAL APPROACH ###
    #####################################
    #####################################

    ##setting up another if statement to handle the situation where orthog == TRUE
    ##in this case, we don't impose the constraint associated with the number of vanishing moments or n.bands
    ##this will result in a potentially dense upsilon and in turn W matrix, but we'll have orthogonal basis functions
    ##the limitation of this approach seems to be that if the basis functions are being built for non-equidistant data, then this is not accounted for in the method

    if(orthog == TRUE){

      ##use the fact that vec(AXB) = vec(C) can be written as (B^T \kronecker A)vec(X) = vec(C) or Mvec(X) = vec(C)
      ##meaning that to find the elements of X we can just use solve()
      mat.A <- Matrix::t(V_j.list[[j]]) %*% V_j.list[[j]]
      mat.B <- Matrix::t(G.tilde_j.list[[j]]) %*% V.tilde_j.list[[j+1]] %*% Matrix::t(V.tilde_j.list[[j+1]]) %*% G.tilde_j.list[[j]]
      mat.M <- Matrix::kronecker(Matrix::t(mat.B), mat.A)

      ##put together C matrix and then make vec(C)
      C.mat <- (-1*Matrix::t(V_j.list[[j]]) %*% V_j.list[[j]] %*% Matrix::t(H.tilde0_j.list[[j]]) %*% V.tilde_j.list[[j+1]] %*% Matrix::t(V.tilde_j.list[[j+1]]) %*% G.tilde_j.list[[j]])
      sol.upsilon_j.list[[j]] <- as(C.mat, "sparseVector")

      ##now that we have the system of equations, let's solve it
      ##solve the system of equations
      el.upsilon_j.solved <- Matrix::solve(mat.M, sol.upsilon_j.list[[j]], tol = 0)

      ##using the row/col info, along with the values of the elements, build a sparse matrix
      upsilon_j.list[[j]] <- Matrix::sparseMatrix(
        i = ((attributes(el.upsilon_j.solved)$i) %% ne_j[[j]]) + 1,
        j = ceiling((attributes(el.upsilon_j.solved)$i+1)/ne_j[[j]]),
        x = attributes(el.upsilon_j.solved)$x,
        dims = c(ne_j[[j]], no_j[[j]]))
    }

    #####################################
    #####################################
    ### CONSTRUCT W.TILDE_j ###
    #####################################
    #####################################

    ##now construct the W.tilde_j matrix
    ##first, obtain the matrix from Eq 46 of Jansen 2016
    eq.46.mat <- rbind(cbind(Matrix::Diagonal(n = ne_j[[j]]), Matrix::Matrix(0, nrow = ne_j[[j]], ncol = no_j[[j]])),
                       cbind(-1*P_j.list[[j]], Matrix::Diagonal(n = no_j[[j]]))) %*%
      rbind(cbind(solve(D_j.list[[j]]), Matrix::Matrix(0, ncol = no_j[[j]], nrow = ne_j[[j]])),
            cbind(Matrix::Matrix(0, ncol = ne_j[[j]], nrow = no_j[[j]]), Matrix::Diagonal(n = no_j[[j]]))) %*%
      rbind(cbind(Matrix::Diagonal(n = ne_j[[j]]), U_j.list[[j]]),
            cbind(Matrix::Matrix(0, nrow = no_j[[j]], ncol = ne_j[[j]]), Matrix::Diagonal(n = no_j[[j]])))

    ##use eq 52 from Jansen 2016 to obtain the W_j matrix - not exactly this eq, but a slightly different version from Jansen's textbook
    W.tilde_j.disjoint[[j]] <- rbind(cbind(Matrix::Diagonal(n = ne_j[[j]]), upsilon_j.list[[j]]),
                                     cbind(Matrix::Matrix(0, ncol = ne_j[[j]], nrow = no_j[[j]]), Matrix::Diagonal(n = no_j[[j]]))) %*%
      eq.46.mat

    ##need to have W.tilde_j in sequential form, not disjoint as it is above
    W.tilde_j[[j]] <- Matrix::Matrix(0, nrow=nrow(W.tilde_j.disjoint[[j]]), ncol=ncol(W.tilde_j.disjoint[[j]]), sparse = TRUE)
    W.tilde_j[[j]][,even.index] <- W.tilde_j.disjoint[[j]][,1:ne_j[[j]]]
    W.tilde_j[[j]][,odd.index] <- W.tilde_j.disjoint[[j]][,(ne_j[[j]]+1):(n_j[[j]])]

    ##now grab H.tilde_j
    H.tilde_j.list[[j]] <- suppressWarnings(Matrix::t(Matrix::Matrix(W.tilde_j[[j]][1:ne_j[[j]],], nrow=ne_j[[j]])))

    ##now expand W.tilde_j to be of the same size for each j
    ident_j <- Matrix::Diagonal(n - n_j[[j]])
    zeroes_j <- Matrix::Matrix(0, nrow = n_j[[j]], ncol = n - n_j[[j]], sparse=TRUE)
    W.tilde_j.expand[[j]] <- rbind(cbind(W.tilde_j[[j]], zeroes_j), cbind(Matrix::t(zeroes_j), ident_j))

    ##get ready for next iteration
    x_j[[j-1]] <- x_j[[j]][even.index]
    n_j[[j-1]] <- length(x_j[[j-1]])
    temp.call[[j]] <- parse(text=paste0("W.tilde_j.expand[[",j,"]]"))
  }

  ##lastly, calculate W.tilde and also the wavelet coefficients from the lifting scheme
  for(i in 2:length(temp.call)){
    if(i == 2){
      temp.expr <- temp.call[[i]]
    } else{
      temp.expr <- paste(temp.expr, temp.call[[i]], sep="%*%")
    }
  }
  W.tilde <- eval(parse(text=temp.expr))

  ##prepare output
  return.objects <- list("x_j"=x_j,
                         "n_j"=n_j,
                         "P_j"=P_j.list,
                         "U_j"=U_j.list,
                         "D_j"=D_j.list,
                         "upsilon_j" = upsilon_j.list,
                         "ne_j"=ne_j,
                         "no_j"=no_j,
                         "H_j" = H_j.list,
                         "G_j" = G_j.list,
                         "H.tilde_j" = H.tilde_j.list,
                         "H.tilde0_j" = H.tilde0_j.list,
                         "G.tilde_j" = G.tilde_j.list,
                         "V_j" = V_j.list,
                         "V.tilde_j" = V.tilde_j.list,
                         "M_j" = M_j.list,
                         "O_j" = O_j.list,
                         "sol.upsilon_j" = sol.upsilon_j.list,
                         "W.tilde_j"= W.tilde_j.expand,
                         "W.tilde"= W.tilde)

  return(return.objects)
}
