import numpy
import scipy
import scipy.linalg

def main():
    eps = 0
    delta = -1

    psi = numpy.zeros(3)
    psi[1] = 1

    dt = 0.1
    t_init = 0.0
    t_final = 50.
    ntimesteps = int((t_final-t_init)/dt + 1)
    times = numpy.linspace(t_init, t_final, ntimesteps)

    omega_c = 0.3
    beta = 1.
    lamda = 0.25
    
    nmodes = 300
    ntraj = 10000


    rho = numpy.zeros((ntimesteps,3,3), dtype=numpy.complex)
    omega_j = numpy.zeros(nmodes)
    
    for i in range(nmodes):

        omega_j[i] = omega_c*numpy.tan((numpy.pi/2.)*(i+0.5)/nmodes)
   
    c_j = omega_j*numpy.sqrt(2*lamda/nmodes)

    for n in range(ntraj):

        q_1 = numpy.random.normal(0.0, 1/(numpy.sqrt(beta)*omega_c), nmodes)

        Hc1 = numpy.dot(c_j,q_1)
        
        q_2 = numpy.random.normal(0.0, 1/(numpy.sqrt(beta)*omega_c), nmodes)

        Hc2 = numpy.dot(c_j,q_2)

        print 'Processing trajectory', n+1,'...'
#        print 'Hc='Hc
        ham = numpy.array([[0.,      0.,     0.],
                           [0.,  Hc1+eps,  delta],
                           [0.,  delta, Hc2+eps]])

        psi_t = numpy.zeros((len(times), 3), dtype=numpy.complex)
        
        for t,time in enumerate(times):
            if t == 0:
                psi_t[0,:] = psi
            else:
                psi_t[t,:] = numpy.dot(scipy.linalg.expm(-(1j)*ham*time), psi)
        
        for i in range(3):
            for j in range(3):
                rho[:,i,j] += numpy.conjugate(psi_t[:,i])*psi_t[:,j]

    print rho

    rho /= ntraj

    with open('3state_site_population5.dat', 'w') as f:
        for (time,rhot) in zip(times,rho):
            f.write('%0.8f %0.8f %0.8f %0.8f\n'%(time,rhot[0,0].real,rhot[1,1].real,rhot[2,2].real))

    print 'Your simulation has completed successfully.'

if __name__ == '__main__':
    main()
