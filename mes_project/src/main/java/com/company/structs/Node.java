package com.company.structs;

public class Node
{
    public double x,y;
    public double temperature; //dodane do agregacji
    public int numOfNode;

    public Node()
    {
        this.x = 0.0D;
        this.y = 0.0D;
    }

    public Node(double x, double y)
    {
        this.x = x;
        this.y = y;
    }

    public String toString()
    {
        return "Node { x =" +
                this.x +
                ", y =" +
                this.y +
                ", numberOfNode =" +
                this.numOfNode +
                " }";
    }
}
