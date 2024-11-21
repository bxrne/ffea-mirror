#include "rod_math_v9.h"
#include "stdio.h"

int main(){
    // regular energy
    rod::float3 e_i = {1,0,0};
    rod::float3 e_i_equil = {1,0,0};
    rod::float3 e_ip1 = {1,0,0};
    rod::float3 m_ip1 = {0,1,0};
    rod::float3 e_ip1_equil = {1,0,0};
    
    float unperturbed_energy = 0;
    unperturbed_energy += rod::get_stretch_energy(1, e_i, e_i_equil);
    unperturbed_energy += rod::get_stretch_energy(1, e_ip1, e_ip1_equil);

    rod::float3 e_i_perturbed = {1.1,0,0};
    rod::float3 e_ip1_perturbed = {0.9,0,0};
    
    float perturbed_stretch_energy = 0;
    perturbed_stretch_energy += rod::get_stretch_energy(1, e_i_perturbed, e_i_equil);
    perturbed_stretch_energy += rod::get_stretch_energy(1, e_ip1_perturbed, e_ip1_equil);

    rod::float3 new_m_ip1;
    // perturb e_ip1 twist
    float twisted_stretch_energy = 0;
    rod::rodrigues_rotation(m_ip1, e_ip1, 0.006283185, new_m_ip1);
    twisted_stretch_energy += rod::get_stretch_energy(1, e_i, e_i_equil);
    twisted_stretch_energy += rod::get_stretch_energy(1, e_ip1, e_ip1_equil);
    if (unperturbed_energy - twisted_stretch_energy < 0.001 && unperturbed_energy - twisted_stretch_energy > -0.001){
        return 0;
    }
    printf("Unperturbed energy: %f", unperturbed_energy);
    printf("Twisted energy: %f", twisted_stretch_energy);
    return 1;
}
