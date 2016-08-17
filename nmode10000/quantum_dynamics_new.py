import numpy
import scipy
import scipy.linalg

def main():
    eps1 = float(input('eps1='))
    eps2 = float(input('eps2='))
    eps3 = float(input('eps3='))
    delta = float(input('delta='))

    # hamiltonian is a 3x3 matrix in this basis
    hamiltonian = numpy.array([[eps1,  delta, 0.0  ],
                               [delta, eps2,  delta],
                               [0.0,   delta, eps3 ]])

    # psi is a length-3 vector
    psi = numpy.zeros(3)
    psi[0] = float(input("psi1="))
    psi[1] = float(input("psi2="))
    psi[2] = float(input("psi3="))
	
    psinorm = numpy.linalg.norm(psi)
    print "psinorm=", psinorm
    
    if numpy.allclose([0.0],[psinorm])==True:
        psi = psi
    else:
        psi = psi/psinorm

    # example of letting the hamiltonian operate on a wavefunction
    # with numpy's 'dot' function
    hampsi = numpy.dot(hamiltonian, psi)
    print "H*psi ="
    print hampsi

    # trying inner product
    psistar = numpy.conjugate(psi).T
    psihampsi = numpy.dot(psistar, hampsi)
    print "innerprod ="
    print psihampsi

    # example to find the eigenvalues and eigenvectors of hamiltonian
    evals, evecs = numpy.linalg.eigh(hamiltonian)
    # notes: 
    # -- hamiltonian is hermitian, so we can use 'eigh' instead of 'eig'
    # -- evecs is a matrix, whose *columns* are the eigenvectors    
    print "evals ="
    print evals
    print "evecs ="
    print evecs
    print ""

    for e_n, psi_n in zip(evals, evecs.transpose()):
        print "eigenvalue =", e_n
        print "eigenvector ="
        print psi_n
        print ""

    # Now try to calculate the time-dependent quantum dynamics
    # of psi(t) = e^{-i H t} psi(0)
	dt = 0.005
	t_init = 0.0
	t_final = 10.0
	ntimesteps = int( (t_final-t_init)/dt + 1 )
    times = numpy.linspace(t_init, t_final, ntimesteps) # start, stop, number of pts
    # in general, psi(t) will be complex!
    psi_t = numpy.zeros((len(times), 3), dtype=numpy.complex)
    for t,time in enumerate(times):
        # using enumerate function:
        # -- 't' is the integer index (like in a for loop)
        # -- 'time' is the float value of time, i.e. times[t]

        if t == 0:
            # use psi from above to initialize psi(t=0)
            psi_t[0,:] = psi
        else:
            # calculate psi(t) here!
            #psi_t[t,:] = numpy.dot(scipy.linalg.expm(-(1j)*hamiltonian*time), psi)
			psi_t[t,:] = ( psi_t[t-1,:]-(1j)*dt*numpy.dot(hamiltonian, psi_t[t-1,:]) - dt**2*numpy.dot(hamiltonian, numpy.dot(hamiltonian, psi_t[t-1,:]))/2.0 )

    t1 = float(input('Enter a time for rho:'))
    
#    if 0.0 > t1 > 10.0:
#        print "t must be between", t_init, "and", t_final, "."
#        t1=0.0
#        print "Setting t equal to zero."
    # Version 1
#   rho11 = numpy.absolute(psi_t[t1,0]) ** 2
#   rho22 = numpy.absolute(psi_t[t1,1]) ** 2
#   rho33 = numpy.absolute(psi_t[t1,2]) ** 2
#   rho12 = numpy.conjugate(psi_t[t1,0])*psi_t[t1,1]
#    rho13 = numpy.conjugate(psi_t[t1,0])*psi_t[t1,2]
#    rho21 = numpy.conjugate(psi_t[t1,1])*psi_t[t1,0]
#    rho23 = numpy.conjugate(psi_t[t1,1])*psi_t[t1,2]
#    rho31 = numpy.conjugate(psi_t[t1,2])*psi_t[t1,0]
#    rho32 = numpy.conjugate(psi_t[t1,2])*psi_t[t1,1]

#    rho = numpy.array([[rho11, rho12, rho13 ],
#                       [rho21, rho22, rho23 ],
#                       [rho31, rho32, rho33 ]])

    # Version 2
    rho = numpy.zeros((3,3,ntimesteps), dtype=numpy.complex)
    for i in range(3):
        for j in range(3):
            rho[i,j] = numpy.conjugate(psi_t[:,i])*psi_t[:,j]

#    # Version 3
#   rho = numpy.zeros((3,3,ntimesteps), dtype=numpy.complex)
#   for t in range(ntimesteps):
#   rho[:,:,t] = numpy.outer(numpy.conjugate(psi_t[t,:]), psi_t[t,:])


    # Version 4
#    rho = numpy.einsum('ti,tj->ijt', numpy.conjugate(psi_t), psi_t)
    
    rhotrace = numpy.trace(rho[:,:,t1]).real

    print "At t=", t1, ", rho=", rho[:,:,t1]

    if numpy.allclose([1.0],[psinorm])==True:
        print "Looks like your wavefunction is normalized!"
    else:
        print "You didn't normalize your wavefunction, did you?  It's OK, I did it for you :).  See?"
    
    print "Trace of rho=", rhotrace
    
    import matplotlib.pyplot as plt
    plt.plot(times, psi_t[:,0].real, times, psi_t[:,1].real, times, psi_t[:,2].real, linewidth=2.0)
    plt.ylabel('Wavefunction Component')
    plt.xlabel('Time')
    plt.axis([0, 10.0, -1.25, 1.25])
    plt.show()

    #plt.plot(times, rho11, times, rho22, times, rho33, linewidth=2.0)
    plt.plot(times, rho[0,0,:], times, rho[1,1,:], times, rho[2,2,:], linewidth=2.0)
    plt.ylabel('Wavefunction Component Probability')
    plt.xlabel('Time')
    plt.axis([0, 10.0, -0.25, 1.25])
    plt.show()


if __name__ == '__main__':
    main()
