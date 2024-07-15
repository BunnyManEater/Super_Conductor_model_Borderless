#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 
#include <iostream>
#include <vector>
#include <math.h>
#include<cmath>
#include<fstream>
#include <string>

long double CalcRadius(long double i_x, long double i_y, long double j_x, long double j_y)
{
    long double r_2 = (i_x - j_x) * (i_x - j_x) + (i_y - j_y) * (i_y - j_y);
    return sqrt(r_2);
}
bool CheckCollision(std::vector<std::pair<long double,long double>>& vortex_coord,long double dzeta_0, long double new_x, long double new_y)
{
    long double distance = 0.0;
    bool collision_avoided = true;
    for (int i = 0; i < vortex_coord.size(); i++)
    {
        distance = CalcRadius(vortex_coord[i].first, vortex_coord[i].second, new_x, new_y);
        if (distance <= dzeta_0) { collision_avoided = false; }
    }
    return collision_avoided;

}
long double Uij(long double pos_x_i, long double pos_y_i, long double pos_x_j, long double pos_y_j, long double lambda_0)
{
    long double rij = sqrt((pos_x_i - pos_x_j)* (pos_x_i - pos_x_j) + (pos_y_i - pos_y_j)* (pos_y_i - pos_y_j));
    long double arg = static_cast<long double>(rij / lambda_0);
    long double rez = std::cyl_bessel_k(0, arg);
    return rez;
}
long double OneVortexImpact(std::vector<std::pair<long double,long double>>& vortex_coord, int vortex_ind, long double vortex_pos_x, long double vortex_pos_y, long double lambda_0)
{
    long double energy = 0.0;
    for (int i = 0; i < vortex_coord.size(); i++)
    {
        if (i != vortex_ind)
        {
            energy += Uij(vortex_pos_x, vortex_pos_y, vortex_coord[i].first, vortex_coord[i].second, lambda_0);
        }
    }
    return energy;
}
long double Wij(long double E_current, long double E_new)
{
    if (E_current >= E_new)
    {
        return 1.;
    }
    else
    {
        long double beta = 30.;
        return exp(-beta * (E_new - E_current));
    }
}
long double CalcStartConfig(std::vector<std::pair<long double,long double>>& vortex_coord, long double lambda_0) //smth wrong
{
    long double energy = 0;
    for (int i = 0; i < vortex_coord.size(); i++)
    {
        for (int j = i + 1; j < vortex_coord.size(); j++)
        {
            energy += Uij(vortex_coord[i].first, vortex_coord[i].second, vortex_coord[j].first, vortex_coord[j].second, lambda_0);
        }
    }
    return energy;
}


int main()
{

    std::ofstream myfile;
    //myfile.open("supercond.txt");
    

    long double pi = 3.1415926535;
    //условия задачи:
    srand(0); //seed for rand()
    long double rect_length = 5000;
    long double rect_height = 5000;

    long double dzeta_0 = 2;
    long double lambda_0 = 180;
    long double move_step = 18;

    std::vector<long double> energy_v;

    std::vector<std::pair <long double,long double>> vortex_coordinates;
    //задаем начальную конфигурацию
    vortex_coordinates.resize(100);
    for (int i = 0; i < 10; i++)
    {
        for (int k = 0; k < 10; k++)
        {
            vortex_coordinates[i *10 + k] = { 2250. + 10. * i,2250. + 10 * k };
        }
    }
    //рачет энергии начальной конфигурации
    long double total_energy = 0.;
    total_energy = CalcStartConfig(vortex_coordinates, lambda_0);
    bool move_possible = true;

    myfile.open( "start.txt");

    for (int i = 0; i < vortex_coordinates.size(); i++)
    {
        myfile << vortex_coordinates[i].first << "," << vortex_coordinates[i].second << "\n";
    }

    myfile.close();

    for (int n = 0; n < 501; n++)  //количество картинок
    {
        for (int k = 0; k < 10*vortex_coordinates.size(); k++) //количество итераций для создания одной картинки
        {
            //выбор случайного вихря и случайного направления
            int shift_vortex_ind = rand() % vortex_coordinates.size(); 
            long double shift_angle = ((long double)rand() / (RAND_MAX)) * 2 * pi;
            long double shift_x = vortex_coordinates[shift_vortex_ind].first + move_step * cos(shift_angle);
            long double shift_y = vortex_coordinates[shift_vortex_ind].second + move_step * sin(shift_angle);
            //периодическое условие: 
            if (shift_x < 0.) { shift_x = rect_length + shift_x; }
            if (shift_x > rect_length) { shift_x = shift_x - rect_length; }
            if (shift_y < 0.) { shift_y = rect_height + shift_y; }
            if (shift_y > rect_height) { shift_y = shift_y - rect_height; }

            //проверка: новая координата находится не в дзета радиусе другого вихря:

            if (CheckCollision(vortex_coordinates, dzeta_0, shift_x, shift_y) == false) {
                move_possible = false;
            }
            if (move_possible == true)
            {
                //расчет энергии новой конфигурации
                long double old_vortex_impact = OneVortexImpact(vortex_coordinates, shift_vortex_ind, vortex_coordinates[shift_vortex_ind].first, vortex_coordinates[shift_vortex_ind].second, lambda_0);
                long double new_vortex_impact = OneVortexImpact(vortex_coordinates, shift_vortex_ind, shift_x, shift_y, lambda_0);
                //расчет вероятности движения 
                long double wij = Wij(old_vortex_impact, new_vortex_impact);
                //Генерация случайного чила (0,1)
                long double R = (double)rand() / (RAND_MAX);
                //Проверка на R<W
                if (R >= wij)
                {
                    move_possible = false;
                }
                if (move_possible == true)
                {
                    vortex_coordinates[shift_vortex_ind].first = shift_x;
                    vortex_coordinates[shift_vortex_ind].second = shift_y;
                    total_energy = total_energy - old_vortex_impact + new_vortex_impact;
                }
            }
            shift_x = 0.;
            shift_y = 0.;
            move_possible = true;
        }
        
        if ( n <= 10 || n == 50 || n == 400 || n == 500 || n == 800 || n== 1000 || n== 2000 || n==3000 || n==4000 || n==4999 ) 
        {
            myfile.open("k" + std::to_string(n) + ".txt");

            for (int i = 0; i < vortex_coordinates.size(); i++)
            {
                myfile << vortex_coordinates[i].first << "," << vortex_coordinates[i].second << "\n";
            }
        }
       
        
        myfile.close();
        energy_v.push_back(total_energy);
    }
    
    for (int i = 0; i < vortex_coordinates.size(); i++)
    {
        std::cout << " Vortex " << i << " x: " << vortex_coordinates[i].first << " y: " << vortex_coordinates[i].second << std::endl;
    }

    long double analyt_energy = CalcStartConfig(vortex_coordinates, lambda_0);

    std::cout << "--------------" << std::endl;
    std::cout << "Iteration energy calc: " << total_energy << "  ||  Analyt_energy: " << analyt_energy << std::endl;

    myfile.open("energy.txt");

    myfile << analyt_energy << "," << total_energy << "\n";
    long double dif = analyt_energy - total_energy;
    //typeid(dif).name();
    myfile.close();
    
    myfile.open("energy_x_y.txt");
    for (int n = 0; n < energy_v.size(); n++)
    {
        myfile << n << "," << energy_v[n] << "\n";
    }
    
    myfile.close();

    return 0;
}

