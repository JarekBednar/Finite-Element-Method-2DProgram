package com.company.structs;

import java.util.Arrays;

public class Jakobian
{
    double[][] jakobian = new double[2][2]; //na jakobian
    double[][] invertedJakobian = new double[2][2]; //na odwrocony jakobian

    double detI; //wyznacznik jakobianu

    public Jakobian(Element elem, Element4_2D elem4_2D, int numOfPnt)
    {
        Node[] nodes = elem.nodes;

        //wyliczenie jakobianu ze wzorow: dx/dKsi = dN1/dKsi*x1 + ... + dN4/dKsi*x4 i adekwatnie dla dy/dEta
        this.jakobian[0][0] = nodes[0].x * elem4_2D.dN_dKsi[numOfPnt][0] + nodes[1].x * elem4_2D.dN_dKsi[numOfPnt][1] + nodes[2].x * elem4_2D.dN_dKsi[numOfPnt][2] + nodes[3].x * elem4_2D.dN_dKsi[numOfPnt][3];
        this.jakobian[0][1] = nodes[0].y * elem4_2D.dN_dKsi[numOfPnt][0] + nodes[1].y * elem4_2D.dN_dKsi[numOfPnt][1] + nodes[2].y * elem4_2D.dN_dKsi[numOfPnt][2] + nodes[3].y * elem4_2D.dN_dKsi[numOfPnt][3];
        this.jakobian[1][0] = nodes[0].x * elem4_2D.dN_dEta[numOfPnt][0] + nodes[1].x * elem4_2D.dN_dEta[numOfPnt][1] + nodes[2].x * elem4_2D.dN_dEta[numOfPnt][2] + nodes[3].x * elem4_2D.dN_dEta[numOfPnt][3];
        this.jakobian[1][1] = nodes[0].y * elem4_2D.dN_dEta[numOfPnt][0] + nodes[1].y * elem4_2D.dN_dEta[numOfPnt][1] + nodes[2].y * elem4_2D.dN_dEta[numOfPnt][2] + nodes[3].y * elem4_2D.dN_dEta[numOfPnt][3];
        //odwrocony jakobian, zapisujemy zamienione z jakobianu zeby sie nie pomieszalo
        this.invertedJakobian[0][0] = this.jakobian[1][1];
        this.invertedJakobian[0][1] = -1.0D * this.jakobian[0][1];
        this.invertedJakobian[1][0] = -1.0D * this.jakobian[1][0];
        this.invertedJakobian[1][1] = this.jakobian[0][0];
        //policzenie wyznacznika
        this.detI = this.jakobian[0][0] * this.jakobian[1][1] - this.jakobian[0][1] * this.jakobian[1][0];
    }

    //druga wersja z zajec, z podawana siatka
    //pochodne z Element4_2D są wykorzystane
    public Jakobian(int numOfElement, int numOfPnt, Element4_2D elem4_2D, Grid grid)
    {
        Node[] nodes = grid.elements[numOfElement].nodes;
        this.jakobian[0][0] = nodes[0].x * elem4_2D.dN_dKsi[numOfPnt][0] + nodes[1].x * elem4_2D.dN_dKsi[numOfPnt][1] + nodes[2].x * elem4_2D.dN_dKsi[numOfPnt][2] + nodes[3].x * elem4_2D.dN_dKsi[numOfPnt][3];
        this.jakobian[0][1] = nodes[0].y * elem4_2D.dN_dKsi[numOfPnt][0] + nodes[1].y * elem4_2D.dN_dKsi[numOfPnt][1] + nodes[2].y * elem4_2D.dN_dKsi[numOfPnt][2] + nodes[3].y * elem4_2D.dN_dKsi[numOfPnt][3];
        this.jakobian[1][0] = nodes[0].x * elem4_2D.dN_dEta[numOfPnt][0] + nodes[1].x * elem4_2D.dN_dEta[numOfPnt][1] + nodes[2].x * elem4_2D.dN_dEta[numOfPnt][2] + nodes[3].x * elem4_2D.dN_dEta[numOfPnt][3];
        this.jakobian[1][1] = nodes[0].y * elem4_2D.dN_dEta[numOfPnt][0] + nodes[1].y * elem4_2D.dN_dEta[numOfPnt][1] + nodes[2].y * elem4_2D.dN_dEta[numOfPnt][2] + nodes[3].y * elem4_2D.dN_dEta[numOfPnt][3];
        //odwrocony jakobian, zapisujemy zamienione z jakobianu zeby sie nie pomieszalo
        this.invertedJakobian[0][0] = this.jakobian[1][1];
        this.invertedJakobian[0][1] = -1.0D * this.jakobian[0][1];
        this.invertedJakobian[1][0] = -1.0D * this.jakobian[1][0];
        this.invertedJakobian[1][1] = this.jakobian[0][0];
        //policzenie wyznacznika
        this.detI = this.jakobian[0][0] * this.jakobian[1][1] - this.jakobian[0][1] * this.jakobian[1][0];
    }

    public String toString()
    {
        String pom = Arrays.toString(this.jakobian[0]);
        String pom2 = "Jakobian { jakobian=[\n" + pom;

        pom2 += ",\n" + Arrays.toString(this.jakobian[1]) +
                "],\njakobianInverse= [\n" + Arrays.toString(this.invertedJakobian[0]) +
                ",\n" + Arrays.toString(this.invertedJakobian[1]) +
                "],\ndet=" + this.detI + "}";

        return pom2;
    }
}
