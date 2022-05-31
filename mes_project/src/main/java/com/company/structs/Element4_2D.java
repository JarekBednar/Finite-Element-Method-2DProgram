package com.company.structs;

import java.util.Arrays;

public class Element4_2D
{
    public double[][] dN_dKsi, dN_dEta; //lokalne, jakby x,y

    public double[][] dN_dx, dN_dy; //macierze do wspolrzednych globalnych
    public double[][] N;

    public int numOfPts; //liczba punktow

    //jako argument liczba punktow
    public Element4_2D(int numOfPts)
    {
        Gauss gauss = new Gauss(numOfPts); //wykorzystanie klasy Gauss
        this.numOfPts = numOfPts; //ilosc punktow
        this.dN_dEta = new double[numOfPts * numOfPts][4]; //wielkosc pierwsza jest podnoszona do kwadratu
        this.dN_dKsi = new double[numOfPts * numOfPts][4]; //wielkosc pierwsza jest podnoszona do kwadratu

        this.dN_dx = null;
        this.dN_dy = null;

        int i;
        int j;

        //uzupelnienie wartosci Ksi w tablicy
        for(i = 0; i < numOfPts; ++i)
        {
            for(j = 0; j < numOfPts; ++j)
            {
                this.dN_dKsi[i * numOfPts + j][0] = -1.0D * (1.0D - gauss.pts[i]) / 4.0D;
                this.dN_dKsi[i * numOfPts + j][1] = (1.0D - gauss.pts[i]) / 4.0D;
                this.dN_dKsi[i * numOfPts + j][2] = -1.0D * (1.0D + gauss.pts[i]) / 4.0D;
                this.dN_dKsi[i * numOfPts + j][3] = (1.0D + gauss.pts[i]) / 4.0D;
            }
        }

        //uzupelnienie wartosci Eta w tablicy
        for(i = 0; i < numOfPts; ++i)
        {
            for(j = 0; j < numOfPts; ++j)
            {
                this.dN_dEta[i * numOfPts + j][0] = -1.0D * (1.0D - gauss.pts[j]) / 4.0D;
                this.dN_dEta[i * numOfPts + j][1] = -1.0D * (1.0D + gauss.pts[j]) / 4.0D;
                this.dN_dEta[i * numOfPts + j][2] = (1.0D - gauss.pts[j]) / 4.0D;
                this.dN_dEta[i * numOfPts + j][3] = (1.0D + gauss.pts[j]) / 4.0D;
            }
        }

        //dla 2-punktowego wyliczenie wartosci funkcji kształtu
        // elementu czterowezlowego w pkt calkowania
        //wzor byl podany 1/4(1-Ksi)(1-Eta)

        if (numOfPts == 2)
        {
            this.N = new double[4][4]; //4,4 bo 2^2,4
            this.N[0][0] = 0.25D * (1.0D - gauss.pts[0]) * (1.0D - gauss.pts[0]);
            this.N[0][1] = 0.25D * (1.0D + gauss.pts[0]) * (1.0D - gauss.pts[0]);
            this.N[0][2] = 0.25D * (1.0D - gauss.pts[0]) * (1.0D + gauss.pts[0]);
            this.N[0][3] = 0.25D * (1.0D + gauss.pts[0]) * (1.0D + gauss.pts[0]);
            this.N[1][0] = 0.25D * (1.0D - gauss.pts[1]) * (1.0D - gauss.pts[0]);
            this.N[1][1] = 0.25D * (1.0D + gauss.pts[1]) * (1.0D - gauss.pts[0]);
            this.N[1][2] = 0.25D * (1.0D - gauss.pts[1]) * (1.0D + gauss.pts[0]);
            this.N[1][3] = 0.25D * (1.0D + gauss.pts[1]) * (1.0D + gauss.pts[0]);
            this.N[2][0] = 0.25D * (1.0D - gauss.pts[0]) * (1.0D - gauss.pts[1]);
            this.N[2][1] = 0.25D * (1.0D + gauss.pts[0]) * (1.0D - gauss.pts[1]);
            this.N[2][2] = 0.25D * (1.0D - gauss.pts[0]) * (1.0D + gauss.pts[1]);
            this.N[2][3] = 0.25D * (1.0D + gauss.pts[0]) * (1.0D + gauss.pts[1]);
            this.N[3][0] = 0.25D * (1.0D - gauss.pts[1]) * (1.0D - gauss.pts[1]);
            this.N[3][1] = 0.25D * (1.0D + gauss.pts[1]) * (1.0D - gauss.pts[1]);
            this.N[3][2] = 0.25D * (1.0D - gauss.pts[1]) * (1.0D + gauss.pts[1]);
            this.N[3][3] = 0.25D * (1.0D + gauss.pts[1]) * (1.0D + gauss.pts[1]);
        }

        //to samo, tylko dla 3-punktowego wyliczenie wartosci funkcji kształtu
        // elementu czterowezlowego w pkt calkowania
        if (numOfPts == 3)
        {
            this.N = new double[9][4]; //9,4 bo 3^2,4
            this.N[0][0] = 0.25D * (1.0D - gauss.pts[0]) * (1.0D - gauss.pts[0]);
            this.N[0][1] = 0.25D * (1.0D + gauss.pts[0]) * (1.0D - gauss.pts[0]);
            this.N[0][2] = 0.25D * (1.0D - gauss.pts[0]) * (1.0D + gauss.pts[0]);
            this.N[0][3] = 0.25D * (1.0D + gauss.pts[0]) * (1.0D + gauss.pts[0]);
            this.N[1][0] = 0.25D * (1.0D - gauss.pts[1]) * (1.0D - gauss.pts[0]);
            this.N[1][1] = 0.25D * (1.0D + gauss.pts[1]) * (1.0D - gauss.pts[0]);
            this.N[1][2] = 0.25D * (1.0D - gauss.pts[1]) * (1.0D + gauss.pts[0]);
            this.N[1][3] = 0.25D * (1.0D + gauss.pts[1]) * (1.0D + gauss.pts[0]);
            this.N[2][0] = 0.25D * (1.0D - gauss.pts[2]) * (1.0D - gauss.pts[0]);
            this.N[2][1] = 0.25D * (1.0D + gauss.pts[2]) * (1.0D - gauss.pts[0]);
            this.N[2][2] = 0.25D * (1.0D - gauss.pts[2]) * (1.0D + gauss.pts[0]);
            this.N[2][3] = 0.25D * (1.0D + gauss.pts[2]) * (1.0D + gauss.pts[0]);
            this.N[3][0] = 0.25D * (1.0D - gauss.pts[0]) * (1.0D - gauss.pts[1]);
            this.N[3][1] = 0.25D * (1.0D + gauss.pts[0]) * (1.0D - gauss.pts[1]);
            this.N[3][2] = 0.25D * (1.0D - gauss.pts[0]) * (1.0D + gauss.pts[1]);
            this.N[3][3] = 0.25D * (1.0D + gauss.pts[0]) * (1.0D + gauss.pts[1]);
            this.N[4][0] = 0.25D * (1.0D - gauss.pts[1]) * (1.0D - gauss.pts[1]);
            this.N[4][1] = 0.25D * (1.0D + gauss.pts[1]) * (1.0D - gauss.pts[1]);
            this.N[4][2] = 0.25D * (1.0D - gauss.pts[1]) * (1.0D + gauss.pts[1]);
            this.N[4][3] = 0.25D * (1.0D + gauss.pts[1]) * (1.0D + gauss.pts[1]);
            this.N[5][0] = 0.25D * (1.0D - gauss.pts[2]) * (1.0D - gauss.pts[1]);
            this.N[5][1] = 0.25D * (1.0D + gauss.pts[2]) * (1.0D - gauss.pts[1]);
            this.N[5][2] = 0.25D * (1.0D - gauss.pts[2]) * (1.0D + gauss.pts[1]);
            this.N[5][3] = 0.25D * (1.0D + gauss.pts[2]) * (1.0D + gauss.pts[1]);
            this.N[6][0] = 0.25D * (1.0D - gauss.pts[0]) * (1.0D - gauss.pts[2]);
            this.N[6][1] = 0.25D * (1.0D + gauss.pts[0]) * (1.0D - gauss.pts[2]);
            this.N[6][2] = 0.25D * (1.0D - gauss.pts[0]) * (1.0D + gauss.pts[2]);
            this.N[6][3] = 0.25D * (1.0D + gauss.pts[0]) * (1.0D + gauss.pts[2]);
            this.N[7][0] = 0.25D * (1.0D - gauss.pts[1]) * (1.0D - gauss.pts[2]);
            this.N[7][1] = 0.25D * (1.0D + gauss.pts[1]) * (1.0D - gauss.pts[2]);
            this.N[7][2] = 0.25D * (1.0D - gauss.pts[1]) * (1.0D + gauss.pts[2]);
            this.N[7][3] = 0.25D * (1.0D + gauss.pts[1]) * (1.0D + gauss.pts[2]);
            this.N[8][0] = 0.25D * (1.0D - gauss.pts[2]) * (1.0D - gauss.pts[2]);
            this.N[8][1] = 0.25D * (1.0D + gauss.pts[2]) * (1.0D - gauss.pts[2]);
            this.N[8][2] = 0.25D * (1.0D - gauss.pts[2]) * (1.0D + gauss.pts[2]);
            this.N[8][3] = 0.25D * (1.0D + gauss.pts[2]) * (1.0D + gauss.pts[2]);
        }

    }


//pochodne potrzebne do wyliczenia macierzy H
    public void count_Derivatives_XY(Jakobian jakob)
    {
        this.dN_dx = new double[this.numOfPts * this.numOfPts][4];

        this.dN_dy = new double[this.numOfPts * this.numOfPts][4];

        //mnożona jest macierz przez wektor
        double[] var = jakob.jakobian[0];
        var[0] *= 1.0D / jakob.detI;
        var = jakob.jakobian[0];
        var[1] *= 1.0D / jakob.detI;
        var = jakob.jakobian[1];
        var[0] *= 1.0D / jakob.detI;
        var = jakob.jakobian[1];
        var[1] *= 1.0D / jakob.detI;

        //dalsze obliczenia z punktami całkowania
        //pc - punkt calkowania
        for(int pc = 0; pc < this.numOfPts * this.numOfPts; ++pc)
        {
            for(int i = 0; i < 4; ++i) {
                //wartosc jakobianu * pc z Ksi + wartosc jakobianu * pc z Eta
                this.dN_dx[pc][i] = jakob.jakobian[0][0] * this.dN_dKsi[pc][i] + jakob.jakobian[0][1] * this.dN_dEta[pc][i];
                this.dN_dy[pc][i] = jakob.jakobian[1][0] * this.dN_dKsi[pc][i] + jakob.jakobian[1][1] * this.dN_dEta[pc][i];
            }
        }

    }

    public String toString()
    {
        String as_String = "Element4_2D { dN_dKsi={\n";

        int i;

        for(i = 0; i < this.numOfPts * this.numOfPts; ++i)
        {
            as_String = as_String + "\t" + Arrays.toString(this.dN_dKsi[i]) + "\n";
        }

        as_String = as_String + " }, dN_dEta = {\n";

        for(i = 0; i < this.numOfPts * this.numOfPts; ++i)
        {
            as_String = as_String + "\t" + Arrays.toString(this.dN_dEta[i]) + "\n";
        }

        if (this.dN_dx != null)
        {
            as_String = as_String + "}, dN_dx = {\n";

            for(i = 0; i < this.numOfPts * this.numOfPts; ++i)
            {
                as_String = as_String + "\t" + Arrays.toString(this.dN_dx[i]) + "\n";
            }
        }

        if (this.dN_dy != null)
        {
            as_String = as_String + "}, dN_dy = {\n";

            for(i = 0; i < this.numOfPts * this.numOfPts; ++i)
            {
                as_String = as_String + "\t" + Arrays.toString(this.dN_dy[i]) + "\n";
            }
        }

        as_String = as_String + " }}";
        return as_String;
    }
}
