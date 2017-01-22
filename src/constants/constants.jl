#! /usr/bin/env julia
"""

        This type contains a description of properties of the LBM simulation
        """
immutable LBM_Constants
    
    U::Float64   # Maximal Velocity
    Re::Float64  # Reynolds-Number
    char_L::Float64 # Characteristic length
    nu::Float64  # Viscocity (LBM)
    nu_F::Float64 # Viscosity (fluid)
    rho_F::Float64
    tau::Float64 
    phys_x::Float64 # Physical dimension of x
    phys_y::Float64
    
    LBM_Constants(_U, _char_L, nu_F, rho_F, _phyx, _phyy, _tau) =
        (
            _nu = (2. * _tau - 1.) / 6.;
            _Re = _char_L * _U / _nu;
            new(_U, _Re, _char_L, _nu, nu_F, rho_F,  _tau,
                _phyx, _phyy)
        )
end

=======
export LBM_Constants
