#include "rod_structure.h"
#include <iostream>

int main(){
    float nodes[9] = {0,0,0, 1,0,0, 2,0,0};
    float current_x[3] = {nodes[0], nodes[1], nodes[2]};
    float current_xp1[3] = {nodes[3], nodes[4], nodes[5]};
    float current_xp2[3] = {nodes[6], nodes[7], nodes[8]};
    float e_i_equil[3];
    float e_ip1_equil[3];
    rod::get_p_i(current_x, current_xp1, e_i_equil);
    rod::get_p_i(current_xp1, current_xp2, e_ip1_equil);
    
    float m_i[3] = {0,1,0};
    float m_ip1[3] = {0,1,0};
    
    nodes[5] += 0.5;
    
    float new_x[3] = {nodes[0], nodes[1], nodes[2]};
    float new_xp1[3] = {nodes[3], nodes[4], nodes[5]};
    float new_xp2[3] = {nodes[6], nodes[7], nodes[8]};
    float e_i[3];
    float e_ip1[3];
    rod::get_p_i(new_x, new_xp1, e_i);
    rod::get_p_i(new_xp1, new_xp2, e_ip1);
    
    float m_i_prime[3];
    float m_ip1_prime[3];
    
    float t_i[3];
    float t_ip1[3];
    float t_i_equil[3];
    float t_ip1_equil[3];    
    rod::normalize(e_i, t_i);
    rod::normalize(e_i_equil, t_i_equil);
    rod::normalize(e_ip1_equil, t_ip1_equil);
    rod::normalize(e_ip1, t_ip1);
    
    rod::update_m1_matrix( m_i, e_i_equil, e_i, m_i_prime );
    rod::update_m1_matrix( m_ip1, e_ip1_equil, e_ip1, m_ip1_prime );
       
    float twist = rod::get_twist_energy(1, m_ip1_prime, m_i_prime, m_ip1, m_i, e_i, e_ip1, e_i_equil, e_ip1_equil);
    
    if (twist<0.001){
        return 0;
    }
    std::cout << twist << "\n";
    return 1;
    // test: doing this material frame update (check on wolfram?) SHOULD NOT MAKE A TWIST ENERGY    
}
