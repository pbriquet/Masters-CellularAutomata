using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OpenTK.Graphics.OpenGL;
using Article1.PiLinearAlgebra;
using System.Drawing;

namespace Article1
{
    static class Visualization
    {

        // 1. Build Visualization Functions

        // 1.1 Points Visualization in OpenGL
        // 1.2 Grid Square Visualization in OpenGL
        // 1.3 Growth Polygons Visualization in OpenGL
        // 1.4 State Cell Visualization in OpenGL
        // 1.5 Descentralization from Polygon Visualization in OpenGL

        // 2. Refresh Functions
 
        public static void PointsVisualization(Vector2 [,] Mesh, int mxp, int myp)
        {
            GL.Color3(Color.Black);

            for (int i = 0; i < myp; i++)
            {
                for (int j = 0; j < mxp; j++)
                {
                    GL.Begin(BeginMode.Points);
                    GL.Vertex2(Mesh[i, j].x, Mesh[i, j].y);
                    GL.End();
                }
            }
        }

        public static void SquaresVisualization(Vector2[,] Mesh, int mxp, int myp, double dxGL)
        {
            GL.Color3(Color.Black);

            for (int i = 0; i < myp; i++)
            {
                for (int j = 0; j < mxp; j++)
                {
                    GL.Begin(BeginMode.LineLoop);
                    GL.Vertex2(Mesh[i, j].x - dxGL / 2.0, Mesh[i, j].y - dxGL / 2.0);
                    GL.Vertex2(Mesh[i, j].x + dxGL / 2.0, Mesh[i, j].y - dxGL / 2.0);
                    GL.Vertex2(Mesh[i, j].x + dxGL / 2.0, Mesh[i, j].y + dxGL / 2.0);
                    GL.Vertex2(Mesh[i, j].x - dxGL / 2.0, Mesh[i, j].y + dxGL / 2.0);
                    GL.End();
                }
            }
        }
        public static void CAPolygonVisualization(Vector2[,][] Polygons, int[,] state, int mxp, int myp)
        {
            for (int i = 0; i < myp; i++)
            {
                for (int j = 0; j < mxp; j++)
                {
                    if (state[i, j] != 0 && Polygons[i, j] != null)
                    {

                        GL.Color3(Color.Green);

                        GL.Begin(BeginMode.LineLoop);

                        for (int k = 0; k < 4; k++)
                        {
                            GL.Vertex2(Polygons[i, j][k].x, Polygons[i, j][k].y);
                        }

                        GL.End();
                    }
                }
            }
        }
        public static void CAstateVisualization(Vector2[,] Mesh, int[,] state, int mxp, int myp, double dxGL)
        {
            for (int i = 0; i < myp; i++)
            {
                for (int j = 0; j < mxp; j++)
                {

                    if (state[i, j] == 1)
                    {
                        GL.Color3(Color.Blue);

                        GL.Begin(BeginMode.Quads);
                        GL.Vertex2(Mesh[i, j].x - 0.5 * dxGL, Mesh[i, j].y - 0.5 * dxGL);
                        GL.Vertex2(Mesh[i, j].x + 0.5 * dxGL, Mesh[i, j].y - 0.5 * dxGL);
                        GL.Vertex2(Mesh[i, j].x + 0.5 * dxGL, Mesh[i, j].y + 0.5 * dxGL);
                        GL.Vertex2(Mesh[i, j].x - 0.5 * dxGL, Mesh[i, j].y + 0.5 * dxGL);
                        GL.End();

                    }
                }
            }
        }

        public static void dCVisualization(Vector2[,] Mesh, Vector2[,] dC, int[,] state, int mxp, int myp)
        {
            for (int i = 0; i < myp; i++)
            {
                for (int j = 0; j < mxp; j++)
                {

                    if (state[i, j] > 0)
                    {
                        GL.Color3(Color.Black);

                        GL.Begin(BeginMode.Lines);
                        GL.Vertex2(Mesh[i, j].x, Mesh[i, j].y);
                        GL.Vertex2(Mesh[i, j].x + dC[i, j].x, Mesh[i, j].y + dC[i, j].y);
                        GL.End();

                        GL.Color3(Color.Brown);
                        GL.Begin(BeginMode.Points);
                        GL.Vertex2(Mesh[i, j].x + dC[i, j].x, Mesh[i, j].y + dC[i, j].y);
                        GL.End();

                    }
                }
            }
        }




        // 2. Refresh Functions

        public static void RefreshStates(RappazGandinCACell[,] CA, int[,] state, int mx, int my, ref bool started)
        {
            //state = new int[my, mx];

            for (int i = 0; i < my; i++)
            {
                for (int j = 0; j < mx; j++)
                {
                    state[i, j] = CA[i, j].state;
                }
            }
            started = true;
        }

        public static void RefreshPolygons(RappazGandinCACell[,] CA, Vector2[,] Mesh, Vector2[,][] Polygons, int[,] state, int mx, int my, double Lxp, double Lyp, double dxGL)
        {
            Vector2 tmp = new Vector2();

            for (int i = 0; i < my; i++)
            {
                for (int j = 0; j < mx; j++)
                {
                    if (state[i, j] > 0)
                    {
                        for (int k = 0; k < 4; k++)
                        {
                            Polygons[i, j][k] = tmp;

                            tmp = new Vector2(CA[i, j].L[k], CA[i, j].theta);
                            tmp.Rotate(new Angle(k * 90.0, "deg"));
                            tmp = tmp + CA[i, j].dC;
                            tmp = (dxGL / (Lyp / mx)) * tmp;
                            tmp = tmp + Mesh[i, j];

                            Polygons[i, j][k] = tmp;
                        }
                    }
                    else
                    {
                        for (int k = 0; k < 4; k++)
                            Polygons[i, j][k] = new Vector2();
                    }

                }
            }
        }
        public static void RefreshdC(RappazGandinCACell[,] CA, int[,] state, Vector2[,] dC, int mx, int my, double Lxp, double Lyp, double dxGL)
        {
            Vector2 tmp = new Vector2();

            for (int i = 0; i < my; i++)
            {
                for (int j = 0; j < mx; j++)
                {
                    if (state[i, j] > 0)
                    {
                        dC[i, j] = dxGL / (Lxp / mx) * CA[i, j].dC;
                    }

                }
            }
        }
    }
}
