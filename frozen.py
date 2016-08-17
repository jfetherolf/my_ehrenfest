import numpy
import scipy
import scipy.linalg

def main():
    eps = 0
    delta = -1

    psi = numpy.zeros(2)
    psi[1] = 1

    dt = 0.1
    t_init = 0.0
    t_final = 50.
    ntimesteps = int((t_final-t_init)/dt + 1)
    times = numpy.linspace(t_init, t_final, ntimesteps)

    omega_c = 0.3
    beta = 1.
    lamda = 0.25
    
    nmodes = 500
    ntraj = 2000


    rho = numpy.zeros((ntimesteps,2,2), dtype=numpy.complex)
    omega_j = numpy.zeros(nmodes)
    
    for i in range(nmodes):
        omega_j[i] = omega_c*numpy.tan((numpy.pi/(2*nmodes))*(i-0.5))

    c_j = omega_j*numpy.sqrt(lamda/nmodes)

    for n in range(ntraj):

        q_init = numpy.random.normal(0.0, 1/(numpy.sqrt(beta)*omega_c), nmodes)

        Hc = numpy.dot(c_j,q_init)

        print 'Processing trajectory', n+1,'...'

        ham = numpy.array([[Hc+eps, delta],
                           [delta, -Hc-eps]])

        psi_t = numpy.zeros((len(times), 2), dtype=numpy.complex)
        
        for t,time in enumerate(times):
            if t == 0:
                psi_t[0,:] = psi
            else:
                psi_t[t,:] = numpy.dot(scipy.linalg.expm(-(1j)*ham*time), psi)
                #psi_t[t,:] = ( psi_t[t-1,:]-(1j)*dt*numpy.dot(ham, psi_t[t-1,:]) - dt**2*numpy.dot(ham, numpy.dot(ham, psi_t[t-1,:]))/2.0 )
        for i in range(2):
            for j in range(2):
                rho[:,i,j] += numpy.conjugate(psi_t[:,i])*psi_t[:,j]

    rho /= ntraj

    with open('site_population4.dat', 'w') as f:
        for (time,rhot) in zip(times,rho):
            f.write('%0.8f %0.8f %0.8f\n'%(time,rhot[1,1].real,rhot[0,0].real))

    print 'Your simulation has completed successfully.'

if __name__ == '__main__':
    main()
