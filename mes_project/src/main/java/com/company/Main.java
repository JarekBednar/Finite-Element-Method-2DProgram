package com.company;

import com.company.structs.Grid;
import java.util.Arrays;

public class Main
{
    public static void main(String[] args)
    {
        double initial_temp = 100.0D; //temperatura poczatkowa
        double simulation_time = 500.0D; //czas koncowy
        double sim_step_time = 50.0D; //delta tau
        double ambient_temp = 1200.0D; //temperatura otoczenia, czy jest dobra, bo zostala podana dla przykladu
        double alpha = 300.0D; //alfa, wspolczynnik konwekcyjnej wymiany ciepla
        double H = 0.1D; //wysokosc
        double B = 0.1D; //szerokosc
        int nH = 4;     //wezly na wysokosc
        int nB = 4;     //wezly na szerokosc
        double spec_heat = 700.0D; //cieplo wlasciwe
        double k = 25.0D; //wspolczynnik przewodzenia
        double density = 7800.0D; //gestosc
        int numOfPts = 2; //liczba punktow




        Grid siatka = new Grid(H, B, nH, nB); //utworzenie siatki
        double[] tempT = new double[siatka.nN]; //tablica na temperatury

        //zapełnienie tablicy temp poczatkowa
        for (int i = 0; i < siatka.nN; ++i)
        {
            tempT[i] = initial_temp;
        }
        System.out.println("\nminimalna temperatura w siatce(poczatek) : " + Arrays.stream(tempT).min() + "\nmaksymalna temperatura w siatce(poczatek): " + Arrays.stream(tempT).max());

            //petla po czasie tworzaca i wypisujaca siatki
            for (double i = sim_step_time; i <= simulation_time; i += sim_step_time)
            {


                if (i == simulation_time)
                {
                    System.out.println("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
                    System.out.println("\nIteracja: " + i + "\nOstateczna siatka: ");
                } else if (i == sim_step_time)
                {
                    System.out.println("\nIteracja: " + i + "\nSiatka:");
                } else
                {
                    System.out.println("\nIteracja: " + i + "\nSiatka:");
                }
                siatka = new Grid(H, B, nH, nB);      //siatka jest nadpisywana
                siatka.count_Local_H(k, numOfPts);      //obliczenie macierzy H lokalnie
                siatka.count_Hbc_And_P(alpha, ambient_temp);   //obliczanie macierzy HBC i wektora P
                siatka.count_Local_C(density, spec_heat, numOfPts);     //liczenie macierzy C lokalnie
                siatka.count_Global_H_C_P();          //zsumowanie wartosci w wezlach dla H C i P
                siatka.count_Global_H_and_P_with_C(sim_step_time, tempT); //rozszerzenie o macierz C

                tempT = siatka.count_T();

                System.out.println(Arrays.toString(tempT));

                double min = Arrays.stream(tempT).min().isPresent() ? Arrays.stream(tempT).min().getAsDouble() : 0;
                double max = Arrays.stream(tempT).max().isPresent() ? Arrays.stream(tempT).max().getAsDouble() : 0;
                System.out.println("minimalna temperatura w siatce: " + min + "\nmaksymalna temperatura w siatce: " + max);

                if (i == simulation_time)
                {
                    System.out.println("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
                }
            }


    }
}
