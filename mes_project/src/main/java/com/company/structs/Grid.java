package com.company.structs;

import java.util.Arrays;
import org.ejml.equation.Equation;
import org.ejml.simple.SimpleMatrix;

public class Grid
{
    public double H, B; //wysokosc, szerokosc

    public int nH; //liczba wezlow po wysokosci
    public int nB; //liczba wezlow po szerokosci
    public int nN; //Number of Nodes, liczba wezlow
    public int nE; //Num of Elements

    public double dx;
    public double dy;

    //tablica nodow
    public Node[] nodes;
    public Element[] elements;
    public double[][] global_H;
    public double[][] global_C;
    public double[] global_P;

    public Grid(double h, double b, int nH, int nB)
    {
        H = h;
        B = b;
        this.nH = nH; //liczba wezlow w wysokosci
        this.nB = nB; //liczba wezlow w szerokosci

        //policzenie liczby wezlow i elementow
        this.nN = nH * nB;
        this.nE = (nH - 1) * (nB - 1);

        //szerokosc i wysokosc pojedynczej scianki, do liczenia wspolrzednych wezlow
        this.dx = b / (double)(nB - 1);
        this.dy = h / (double)(nH - 1);

        this.global_H = new double[this.nN][this.nN];
        this.global_P = new double[this.nN];
        this.global_C = new double[this.nN][this.nN];

        //tablica nodow rowna wyliczonej liczbie
        this.nodes = new Node[this.nN];

        int i;
        int j;

        //petla w petli idzie po szerokosci nB, co pozwala "przeskoczyc w odpowiednim momencie
        for(j = 0; j < nH; ++j)
        {
            for(i = 0; i < nB; ++i)
            {
                this.nodes[i * nB + j] = new Node();     //stworzenie noda
                this.nodes[i * nB + j].numOfNode = i * nB + j;   //numer noda, od 0
                this.nodes[i * nB + j].x = (double)j * this.dx; //przypisanie wspolrzednych x i y
                this.nodes[i * nB + j].y = (double)i * this.dy; //liczone mnozac dx i dy
            }
        }

        //tablica elementow rowna wyliczonej liczbie
        this.elements = new Element[this.nE];

        //nH-1 i nB-1 no bo to do liczby elementow
        for(i = 0; i < nH - 1; ++i)
        {
            for(j = 0; j < nB - 1; ++j)
            {
                this.elements[i * (nB - 1) + j] = new Element();    //tworzenie elementu
                this.elements[i * (nB - 1) + j].numOfElem = i * (nB - 1) + j; //numer elementu, od 0
                this.elements[i * (nB - 1) + j].nodes[0] = this.nodes[i * nB + j];
                this.elements[i * (nB - 1) + j].nodes[1] = this.nodes[i * nB + j + 1];
                this.elements[i * (nB - 1) + j].nodes[2] = this.nodes[(i + 1) * nB + j];
                this.elements[i * (nB - 1) + j].nodes[3] = this.nodes[(i + 1) * nB + j + 1];
            }
        }

    }

    //funkcja zliczajaca globalnie H C i P dla wezlow
    //agregacja tych macierzy i wektora
    public void count_Global_H_C_P()
    {
        //rozmiar jest duzy, bo liczba nodow na liczbe nodow
        this.global_H = new double[this.nN][this.nN];
        this.global_C = new double[this.nN][this.nN];
        this.global_P = new double[this.nN];

        for(int i = 0; i < this.nE; ++i)
        {
            for(int x = 0; x < 4; ++x)
            {
                for(int y = 0; y < 4; ++y)
                {
                    //H
                    double[] mac_H = this.global_H[this.elements[i].nodes[x].numOfNode];
                    int elem = this.elements[i].nodes[y].numOfNode;
                    mac_H[elem] += this.elements[i].H[x][y] + this.elements[i].Hbc[x][y];
                    //C
                    double[] mac_C = this.global_C[this.elements[i].nodes[x].numOfNode];
                    int elem2 = this.elements[i].nodes[y].numOfNode;
                    mac_C[elem2] += this.elements[i].C[x][y];
                }
                //P
                double[] vek_P = this.global_P;
                int elem3 = this.elements[i].nodes[x].numOfNode;
                vek_P[elem3] += this.elements[i].P[x];
            }
        }
    }

    //liczenie lokalnie H
    public void count_Local_H(double k, int numOfPts)
    {
        Element4_2D element4_2D = new Element4_2D(numOfPts);

        for(int i = 0; i < this.elements.length; ++i)
        {
            this.elements[i].countH(element4_2D, k);
        }

    }


    public void count_Local_C(double ro, double cp, int numOfPts)
    {
        for(int i = 0; i < this.elements.length; ++i)
        {
            this.elements[i].countC(ro, cp, numOfPts);
        }

    }

//HBC i P na podstawie wzorow, ale bezposrednio bez uzywania BC Boundery Condition(warunek brzegowy)
    public void count_Hbc_And_P(double alpha, double temp_Of_Surroundings)
    {
        Gauss gauss = new Gauss(2);
        double[][] N = new double[2][4];

        int i;
        //po elementach
        for(i = 0; i < this.nE; ++i)
        {
            this.elements[i].Hbc = new double[4][4]; //Hbc jest 4x4
            this.elements[i].P = new double[4];
        }

        double[] pom;
        int x;
        double det;
        int y;
        int z;

        //po l. wezlow na szerokosc
        for(i = 0; i < this.nB - 1; ++i)
        {
            for(x = 0; x < 2; ++x)
            {
                N[x][0] = 0.25D * (1.0D - gauss.pts[x]) * 2.0D;
                N[x][1] = 0.25D * (1.0D + gauss.pts[x]) * 2.0D;
                N[x][2] = 0.25D * (1.0D - gauss.pts[x]) * 0.0D;
                N[x][3] = 0.25D * (1.0D + gauss.pts[x]) * 0.0D;

                det = this.dx / 2.0D;

                for(y = 0; y < 4; ++y)
                {
                    for(z = 0; z < 4; ++z)
                    {
                        pom = this.elements[i].Hbc[y];
                        pom[z] += gauss.weights[x] * alpha * N[x][y] * N[x][z] * det;
                    }

                    pom = this.elements[i].P;
                    pom[y] += alpha * gauss.weights[x] * N[x][y] * temp_Of_Surroundings * det;
                }
            }

        }

        for(i = this.nB - 2; i < this.nE; i += this.nB - 1)
        {
            for(x = 0; x < 2; ++x)
            {
                N[x][0] = 0.0D * (1.0D - gauss.pts[x]);
                N[x][1] = 0.5D * (1.0D - gauss.pts[x]);
                N[x][2] = 0.0D * (1.0D + gauss.pts[x]);
                N[x][3] = 0.5D * (1.0D + gauss.pts[x]);

                det = this.dy / 2.0D;

                for(y = 0; y < 4; ++y)
                {
                    for(z = 0; z < 4; ++z)
                    {
                        pom = this.elements[i].Hbc[y];
                        pom[z] += gauss.weights[x] * alpha * N[x][y] * N[x][z] * det;
                    }

                    pom = this.elements[i].P;
                    pom[y] += alpha * gauss.weights[x] * N[x][y] * temp_Of_Surroundings * det;
                }
            }


        }

        for(i = this.nE - 1; i > this.nE - this.nB; --i)
        {
            for(x = 0; x < 2; ++x)
            {
                N[x][0] = 0.25D * (1.0D - gauss.pts[x]) * 0.0D;
                N[x][1] = 0.25D * (1.0D + gauss.pts[x]) * 0.0D;
                N[x][2] = 0.25D * (1.0D - gauss.pts[x]) * 2.0D;
                N[x][3] = 0.25D * (1.0D + gauss.pts[x]) * 2.0D;

                det = this.dx / 2.0D;

                for(y = 0; y < 4; ++y)
                {
                    for(z = 0; z < 4; ++z)
                    {
                        pom = this.elements[i].Hbc[y];
                        pom[z] += gauss.weights[x] * alpha * N[x][y] * N[x][z] * det;
                    }

                    pom = this.elements[i].P;
                    pom[y] += alpha * gauss.weights[x] * N[x][y] * temp_Of_Surroundings * det;
                }
            }


        }

        for(i = this.nE - (this.nB - 1); i >= 0; i -= this.nB - 1)
        {
            for(x = 0; x < 2; ++x)
            {
                N[x][0] = 0.5D * (1.0D - gauss.pts[x]);
                N[x][1] = 0.0D * (1.0D - gauss.pts[x]);
                N[x][2] = 0.5D * (1.0D + gauss.pts[x]);
                N[x][3] = 0.0D * (1.0D + gauss.pts[x]);

                det = this.dy / 2.0D;

                for(y = 0; y < 4; ++y)
                {
                    for(z = 0; z < 4; ++z)
                    {
                        pom = this.elements[i].Hbc[y];
                        pom[z] += gauss.weights[x] * alpha * N[x][y] * N[x][z] * det;
                    }

                    pom = this.elements[i].P;
                    pom[y] += alpha * gauss.weights[x] * N[x][y] * temp_Of_Surroundings * det;
                }
            }

        }

    }


    public void count_Global_H_and_P_with_C(double dt, double[] T0)
    {
        for(int x = 0; x < this.nN; ++x)
        {
            for(int y = 0; y < this.nN; ++y)
            {
                double[] pom = this.global_H[x];
                pom[y] += this.global_C[x][y] / dt;
                pom = this.global_P;
                pom[x] += this.global_C[x][y] * T0[y] / dt;
            }
        }

    }

    //drukowanie macierzy C
    public void printGlobalHandC()
    {
        String pom = new String();


        for(int i = 0; i < this.nN; ++i)
        {
            pom = pom + Arrays.toString(this.global_H[i]) + "\n";
        }

        pom = pom + "}, globalC = {\n";

        for(int i = 0; i < this.nN; ++i)
        {
            pom = pom + Arrays.toString(this.global_C[i]) + "\n";
        }

        System.out.println(pom);
    }

    //drukowanie wektora
    public void printGlobalP()
    {
        System.out.println(Arrays.toString(this.global_P));
    }

    //liczenie temperatur za pomoca zewnetrznych funkcji
    public double[] count_T()
    {
        double[] T = new double[this.nN];
        SimpleMatrix simple_Matrix1 = new SimpleMatrix(this.global_H);
        SimpleMatrix simple_Matrix2 = new SimpleMatrix(this.nN, 1);

        for(int i = 0; i < this.nN; ++i)
        {
            simple_Matrix2.setRow(i, 0, new double[]{this.global_P[i]});
        }

        Equation equation = new Equation();
        equation.alias(new Object[]{
                simple_Matrix1, "H", simple_Matrix2, "P"});
        equation.process("t=P/H");

        try
        {
            SimpleMatrix solution = (SimpleMatrix)simple_Matrix1.solve(simple_Matrix2);

            for(int i = 0; i < this.nN; ++i)
            {
                T[i] = solution.get(i, 0);
            }
        } catch (Exception var7)
        {
            System.out.println("Blad przy liczeniu temperatur\n");
            var7.printStackTrace();
        }

        return T;
    }

    public String toString()
    {
        System.out.println("Nodes: ");
        Node[] var1 = this.nodes;
        int i = var1.length;

        int var3;
        for(var3 = 0; var3 < i; ++var3)
        {
            Node temp = var1[var3];
            System.out.println("\t" + temp + "\n");
        }

        System.out.println("Elements:   ");
        Element[] var5 = this.elements;
        i = var5.length;

        for(var3 = 0; var3 < i; ++var3)
        {
            Element temp = var5[var3];
            System.out.println("\t" + temp + "\n");
        }

        String temp = new String();

        for(i = 0; i < this.nN; ++i)
        {
            temp = temp + Arrays.toString(this.global_H[i]) + "\n";
        }

        return "Grid { B = " + this.B +
                ", H = " + this.H +
                ", nB = " + this.nB +
                ", nH = " + this.nH +
                ", nE = " + this.nE +
                ", nN = " + this.nN +
                ", dx = " + this.dx +
                ", dy = " + this.dy +
                ", globalH = {\n" + temp + "\n}}";
    }
}
