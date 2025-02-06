include("sde-fcts.jl")

#############################
#
# SETTING UP MODEL PARAMETERS
#
#############################


# this contains the model parameters
#   and running it "empty" assigns
#   default parameter values
# hovering over PARS() should give a summary
p = PARS()

# you can change parameter values as follows
p.K = 20;
p.S = 10;
p.r = zeros(p.S,p.K);

# or you can set them while calling the constructor
p = PARS(K=20,S=10);

#
# TODO: define p by loading in a csv file ...
#




#######################
#
# RUN THE SIMULATION !
#
#######################


# the simulation is ran simply as follows``
out = runSim(p)

# NOTE:
#   if it says "DtLessThanMin"
#   then you probably made some competition coeffs negative
#   and a mutualistic orgy is taking place
#   julia can't handle infinite abundances ... so it aborts




#######################
#
# PLOTTING THE OUTPUT
#
#######################


theme(:orange) # i just love this theme

# you can plot the dynamics of all microbe abundances simultaneously with
plot(out,legend=false,background_color=:transparent,background_color_inside=:transparent)





#######################
#
# SAVING THE PARAMETERS
#
#######################


# this saves the model parameters to a csv file
#       using default file names
savePARS(p)

# you can specify filenames like
savePARS(p,rbN0="rbN0_i")

# or like
savePARS(p,TKS="SOMETHING_COMPLETELY_DIFFERENT",alpha="AlpHa",enviro="🐛🍔🌞")


########################
#
# SAVING THE TIME SERIES
#
########################


# this writes the time series data to a csv file named "abunds.csv"
#   also provides a csv named "times.csv" with time points as a column
saveTS(out)

# filenames can be changed as
saveTS(out,abunds="N",times="t")

# temporal resolution can also be changed like
saveTS(out,pts=10)

# trying for more time points than whats available
#   just defaults to max time points available
saveTS(out,pts=1e100)


###############################
#
# A NOTE ON MATRICIZING VECTORS
#
###############################

# data are written to csv's in "vectorized" form
# that means if something is a matrix (like r or α[k])
# then it is first made into vector where its columns are
# stacked on top of each other, and then this vector is 
# made a column in the associated csv file

# so you might want to take one of these vectorized matrices
# while in R or python or ... and turn it back into a matrix

# below is some code for doing this in julia, where the original
# matrix is assumed to have S rows and K columns... hope this is helpful!


# matricize a vector that was once a S×K matrix
function mat(v,S,K)

    # check if this makes sense
    if length(v) != S*K
        println("DIMENSION MISMATCH ☹️")
    end

    # do the thing
    m = zeros(S,K)
    for k in 1:K
        m[:,k] = v[1+(k-1)*S:k*S]
    end

    return m

end