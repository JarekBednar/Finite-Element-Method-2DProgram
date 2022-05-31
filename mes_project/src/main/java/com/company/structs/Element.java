package com.company.structs;

import java.util.Arrays;

public class Element
{
    //poczatek
    public int[] id;
    public int numOfElem;
    public Node[] nodes;

    //dodane pozniej
    public double[][] H;
    public double[][] C;
    public double[][] Hbc;
    public double[] P;

    public Element()
    {
        this.id = new int[4];
        this.nodes = new Node[4];
        //dodane ostatnio
        this.Hbc = new double[4][4];
        this.P = new double[4];
    }

    public Element(int[] id)
    {
        this.id = id;
    }

    //wyliczanie pojedynczej macierzy H ze wzoru
    public void countH(Element4_2D elem4_2D, double k)
    {
        this.H = new double[4][4]; //macierz H to macierz 4x4
        Gauss gauss = new Gauss(elem4_2D.numOfPts);

        for(int pc = 0; pc < elem4_2D.numOfPts * elem4_2D.numOfPts; ++pc)
        {
            Jakobian jakobian = new Jakobian(this, elem4_2D, pc);
            elem4_2D.count_Derivatives_XY(jakobian);

            for(int i = 0; i < 4; ++i)
            {
                for(int j = 0; j < 4; ++j)
                {
                    double[] suma = this.H[i];
                    //podstawienie pod wzor
                    //zsumowanie macierzy
                    suma[j] += jakobian.detI * k * gauss.weights[pc % elem4_2D.numOfPts] * (elem4_2D.dN_dx[pc][i] * elem4_2D.dN_dx[pc][j] + elem4_2D.dN_dy[pc][i] * elem4_2D.dN_dy[pc][j]);
                }
            }
        }

    }

    public void countC(double ro, double cp, int numOfPts) {
        Gauss gauss = new Gauss(numOfPts);
        Element4_2D elem4_2D = new Element4_2D(numOfPts);
        double[][] N = elem4_2D.N;
        double[][] C = new double[4][4];

        for(int p = 0; p < numOfPts * numOfPts; ++p) {
            Jakobian jakobian = new Jakobian(this, elem4_2D, p);

            for(int x = 0; x < 4; ++x)
            {
                for(int y = 0; y < 4; ++y)
                {
                    //wzor macierzy C
                    C[x][y] += N[p][x] * N[p][y] * ro * cp * jakobian.detI * gauss.weights[p % numOfPts];
                }
            }
        }

        //przypisanie Elementowi
        this.C = C;
    }

    public String toString()
    {
        String pom = "Element { nodes = " + Arrays.toString(this.nodes);
        int i;
        if (this.H != null)
        {
            pom = pom + "}, H = {\n";

            for(i = 0; i < 4; ++i)
            {
                pom = pom + "\t" + Arrays.toString(this.H[i]) + "\n";
            }
        }

        pom = pom + "}, Hbc = {\n";

        for(i = 0; i < 4; ++i)
        {
            pom = pom + Arrays.toString(this.Hbc[i]) + ",\n";
        }

        pom = pom + "P = {\n";
        pom = pom + Arrays.toString(this.P) + ",\n";
        pom = pom + "}, numOfElem = " + this.numOfElem + " }";
        return pom;
    }
}
