package com.company.structs;

import java.util.Arrays;

public class Gauss
{
    public double[] pts; //tablica na punkty
    public double[] weights; //tablica na wagi

    public Gauss(int numOfPts)
    {
        this.pts = new double[numOfPts]; //tablica na punkty
        this.weights = new double[numOfPts]; //za duze, ale na pewno wystarczy

        //jesli sa 2 punkty
        if (numOfPts == 2)
        {
            //odleglosci i wagi, przypisanie odpowiednim punktom odpowiednich wag
            this.pts[0] = -1.0D / Math.sqrt(3.0D);
            this.pts[1] = 1.0D / Math.sqrt(3.0D);
            this.weights[0] = 1.0D;
            this.weights[1] = 1.0D;
        }

        //jesli sa 3 punkty
        if (numOfPts == 3)
        {
            //odleglosci i wagi, przypisanie odpowiednim punktom odpowiednich wag
            this.pts[0] = -Math.sqrt(0.6D);
            this.pts[1] = 0.0D;
            this.pts[2] = Math.sqrt(0.6D);
            //wagi
            this.weights[0] = 0.5555555555555556D;
            this.weights[1] = 0.8888888888888888D;
            this.weights[2] = 0.5555555555555556D;
        }

    }

    //te 4 funkcje licza Gaussa na sztywno z 2 zajec
    //wartosci funkcji 1d i 2d dla zadanej funkcji w wybranym x/x,y
    public double fun1D(double x)
    {
        double fun1D = 5.0D * x * x + 3.0D * x + 6.0D;
        return fun1D;
    }
    public double fun2D(double x, double y)
    {
        double fun2D = 5.0D * x * x * y * y + 3.0D * x * y + 6.0D;
        return fun2D;
    }
    //liczenie calki 1D
    public double countIntegral1D(int numOfPts)
    {
        double integral1D = 0.0D;

        for(int x = 0; x < numOfPts; ++x)
        {
            //zwiekszane o waga danego x * wartosc f w punkcie
            integral1D = integral1D + this.weights[x] * this.fun1D(this.pts[x]);
        }

        return integral1D;
    }
    //liczenie calki 2D
    public double countIntegral2(int numOfPts)
    {
        double integral2D = 0.0D;

        for(int x = 0; x < numOfPts; ++x)
        {
            for(int y = 0; y < numOfPts; ++y)
            {
                //zwiekszanie o waga x * waga y * wartosci f w punkcie po wszystkich
                integral2D = integral2D + this.weights[x] * this.weights[y] * this.fun2D(this.pts[x], this.pts[y]);
            }
        }

        return integral2D;
    }

    public String toString()
    {
        String as_String = Arrays.toString(this.pts);
        return "Gauss { points = " + as_String
                + ", weights = " + Arrays.toString(this.weights) + "}";
    }
}
