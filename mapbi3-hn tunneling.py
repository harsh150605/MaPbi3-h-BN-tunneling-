
#Authors: Harshvardhan motla
#Institution: Delhi university
#Date: July 2025

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, e, m_e, N_A

eV_to_J = 1.602176634e-19 
print("MAPbI₃/h-BN/MoS₂ Heterostructure: Tunneling Probability and H₂ Production Analysis")
print("=" * 85)

# System parameters
BARRIER_HEIGHT_V0 = -1.3      
ELECTRON_ENERGY_E = -3.9     
BARRIER_THICKNESS = 3.3e-10   
EFFECTIVE_MASS = 2.37e-31

# Direct electron flux (independent of photon energy)
ELECTRON_FLUX = 2.7675e17  # electrons/s - direct input parameter
QUANTUM_EFFICIENCY = 0.8

print(f"Input Parameters:")
print(f"├── Barrier height (V₀): {BARRIER_HEIGHT_V0} eV")
print(f"├── Electron energy (E): {ELECTRON_ENERGY_E} eV")
print(f"├── h-BN thickness: {BARRIER_THICKNESS*1e9:.2f} nm")
print(f"├── Effective mass: {EFFECTIVE_MASS/m_e:.2f} mₑ")
print(f"├── Electron flux: {ELECTRON_FLUX:.1e} e⁻/s")
print(f"└── Quantum efficiency: {QUANTUM_EFFICIENCY*100}%")
print()

def calculate_tunneling_probability(V0, E, d, m_eff):
    """
    Calculate quantum tunneling probability through h-BN barrier.
    
    Parameters:
    -----------
    V0 : float
        Barrier height in eV
    E : float  
        Electron energy in eV
    d : float
        Barrier thickness in meters
    m_eff : float
        Effective mass in kg
        
    Returns:
    --------
    tuple : (T, k, barrier_height, electron_energy)
    """
    
    barrier_height = abs(V0 - E)  
    electron_energy = abs(E)      
    
    # Convert to SI units
    barrier_height_J = barrier_height * eV_to_J
    
    # Calculate wave number k = √(2m(V₀-E)/ℏ²)
    k = np.sqrt(2 * m_eff * barrier_height_J) / hbar
    
    # Calculate tunneling probability using exact formula
    # T = [1 + (V₀²sinh²(kd))/(4|E|·|V₀-E|)]⁻¹
    sinh_kd = np.sinh(k * d)
    numerator = (abs(V0)**2 * sinh_kd**2)
    denominator = (4 * electron_energy * barrier_height)
    
    T = 1 / (1 + (numerator / denominator))
    
    return T, k, barrier_height, electron_energy

# Calculate tunneling probability
T, k, barrier_height, electron_energy = calculate_tunneling_probability(
    BARRIER_HEIGHT_V0, ELECTRON_ENERGY_E, BARRIER_THICKNESS, EFFECTIVE_MASS
)

print("1. QUANTUM TUNNELING ANALYSIS")
print("-" * 50)
print(f"Effective barrier height: {barrier_height:.2f} eV")
print(f"Wave number (k): {k:.2e} m⁻¹")
print(f"kd product: {k * BARRIER_THICKNESS:.3f}")
print(f"TUNNELING PROBABILITY: {T:.4f} ({T*100:.2f}%)")
print()

def calculate_electron_transport(electron_flux, quantum_eff, tunneling_prob):
    """
    Calculate electron transport parameters.
    
    Parameters:
    -----------
    electron_flux : float
        Incident electron flux (electrons/s)
    quantum_eff : float
        Internal quantum efficiency (0 to 1)
    tunneling_prob : float
        Tunneling probability (0 to 1)
        
    Returns:
    --------
    dict : Dictionary containing transport parameters
    """
    
    # Effective electron generation considering quantum efficiency
    effective_electrons = electron_flux * quantum_eff
    
    # Electrons that successfully tunnel through h-BN
    tunneling_electrons = effective_electrons * tunneling_prob
    
    # Calculate tunneling current
    tunneling_current = tunneling_electrons * e  # Amperes
    
    return {
        'effective_electrons': effective_electrons,
        'tunneling_electrons': tunneling_electrons,
        'tunneling_current': tunneling_current
    }

# Calculate transport parameters
transport = calculate_electron_transport(ELECTRON_FLUX, QUANTUM_EFFICIENCY, T)

print("2. ELECTRON TRANSPORT ANALYSIS")
print("-" * 50)
print(f"Effective electron flux: {transport['effective_electrons']:.2e} e⁻/s")
print(f"Tunneling electrons: {transport['tunneling_electrons']:.2e} e⁻/s")
print(f"Tunneling current: {transport['tunneling_current']*1e12:.0f} pA")
print()

def calculate_hydrogen_production(tunneling_electron_rate):
    """
    Calculate hydrogen production rate and efficiency.
    
    For water splitting: 2H⁺ + 2e⁻ → H₂
    Each H₂ molecule requires 2 electrons
    
    Parameters:
    -----------
    tunneling_electron_rate : float
        Rate of electrons tunneling through barrier (e⁻/s)
        
    Returns:
    --------
    dict : Dictionary containing hydrogen production parameters
    """
    
    # H₂ production rate (molecules per second)
    h2_molecules_per_sec = tunneling_electron_rate / 2
    
    # Convert to moles
    h2_moles_per_sec = h2_molecules_per_sec / N_A
    h2_moles_per_hour = h2_moles_per_sec * 3600
    h2_moles_per_day = h2_moles_per_hour * 24
    
    # Volume at STP (22.4 L/mol)
    h2_volume_per_day = h2_moles_per_day * 22.4  # L/day
    
    return {
        'h2_molecules_per_sec': h2_molecules_per_sec,
        'h2_moles_per_sec': h2_moles_per_sec,
        'h2_moles_per_hour': h2_moles_per_hour,
        'h2_moles_per_day': h2_moles_per_day,
        'h2_volume_per_day': h2_volume_per_day
    }

# Calculate hydrogen production
hydrogen = calculate_hydrogen_production(transport['tunneling_electrons'])

print("3. HYDROGEN PRODUCTION ANALYSIS")
print("-" * 50)
print(f"H₂ production rate: {hydrogen['h2_molecules_per_sec']:.2e} molecules/s")
print(f"H₂ production rate: {hydrogen['h2_moles_per_day']:.3e} mol/day")
print(f"H₂ volume (STP): {hydrogen['h2_volume_per_day']:.3f} L/day")
print()

def calculate_efficiencies(transport_params, electron_flux):
    """
    Calculate system efficiencies.
    
    Parameters:
    -----------
    transport_params : dict
        Electron transport parameters
    electron_flux : float
        Input electron flux
        
    Returns:
    --------
    dict : Dictionary containing efficiency parameters
    """
    
    # Electron tunneling efficiency
    tunneling_efficiency = transport_params['tunneling_electrons'] / electron_flux
    
    # Current density (assuming 1 cm² active area)
    area = 1e-4  # 1 cm² in m²
    current_density = transport_params['tunneling_current'] / area
    
    return {
        'tunneling_efficiency': tunneling_efficiency,
        'current_density': current_density
    }

# Calculate efficiencies
efficiencies = calculate_efficiencies(transport, ELECTRON_FLUX)

print("4. SYSTEM EFFICIENCY ANALYSIS")
print("-" * 50)
print(f"Electron tunneling efficiency: {efficiencies['tunneling_efficiency']*100:.2f}%")
print(f"Current density (1 cm²): {efficiencies['current_density']:.0f} A/m²")
print()

# Results summary
print("5. SUMMARY OF KEY RESULTS")
print("=" * 85)
print(f"System: MAPbI₃/h-BN/MoS₂ van der Waals heterostructure")
print(f"h-BN barrier: {BARRIER_THICKNESS*1e9:.2f} nm thick, {barrier_height:.1f} eV barrier height")
print()
print(f"Key Performance Metrics:")
print(f"├── Tunneling probability: {T:.3f} ({T*100:.1f}%)")
print(f"├── Tunneling current: {transport['tunneling_current']*1e12:.0f} pA")
print(f"├── H₂ production rate: {hydrogen['h2_molecules_per_sec']:.2e} molecules/s")
print(f"├── H₂ production rate: {hydrogen['h2_moles_per_day']*1000:.1f} mmol/day")
print(f"├── H₂ volume (STP): {hydrogen['h2_volume_per_day']*1000:.0f} mL/day")
print(f"└── Tunneling efficiency: {efficiencies['tunneling_efficiency']*100:.1f}%")
print("=" * 85)

def parametric_analysis():
    """
    Perform parametric analysis for sensitivity study.
    """
    
    # Thickness dependence
    thickness_range = np.linspace(0.1e-10, 1.0e-10, 50)  # 0.01 to 0.1 nm
    T_vs_thickness = []
    
    for d in thickness_range:
        T_param, _, _, _ = calculate_tunneling_probability(
            BARRIER_HEIGHT_V0, ELECTRON_ENERGY_E, d, EFFECTIVE_MASS
        )
        T_vs_thickness.append(T_param)
    
    # Create figure for publication
    plt.figure(figsize=(10, 6))
    plt.plot(thickness_range*1e9, np.array(T_vs_thickness)*100, 'b-', linewidth=2)
    plt.axvline(BARRIER_THICKNESS*1e9, color='red', linestyle='--', 
                label=f'Current design: {BARRIER_THICKNESS*1e9:.2f} nm', linewidth=2)
    plt.xlabel('h-BN Thickness (nm)', fontsize=12)
    plt.ylabel('Tunneling Probability (%)', fontsize=12)
    plt.title('Tunneling Probability vs h-BN Barrier Thickness', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    plt.tight_layout()
    
    # Save figure
    plt.savefig('C:\\Users\\my\\tunneling_focused_analysis.png', dpi=300, bbox_inches='tight')
    print(f"\nParametric analysis plot saved as 'tunneling_focused_analysis.png'")
    
    return thickness_range, T_vs_thickness

# Run parametric analysis
thickness_data, tunneling_data = parametric_analysis()

print(f"\nScript execution completed successfully.")
print(f"Focus: Tunneling probability and hydrogen production analysis.")
print(f"Results suitable for peer-reviewed publication.")
