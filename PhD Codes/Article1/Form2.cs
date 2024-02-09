using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Threading;
using OpenTK.Graphics.OpenGL;
using Article1.PiLinearAlgebra;
using System.Reflection;

namespace Article1
{
    public partial class Form2 : Form
    {

        bool GLloaded = false;
        bool started = false;
        bool pause = false;
        bool buildingVisual = false;
        Vector2[,] Mesh;
        Vector2[,][] Polygons;
        Vector2[,] dC;

        int[,] state;

        int mxp;
        int myp;

        double Lxp;
        double Lyp;

        double LxVisual;
        double LyVisual;

        double dxGL = new double();

        Model.Al MyAlloy = new Model.Al();

        public Form2()
        {
            InitializeComponent();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            if (button1.Text == "Run")
            {
                button1.Text = "Pause";
                RunParameter v = new RunParameter(tmaxtxtbox.Text, dtmaxtxtbox.Text, rtitxtbox.Text, mxtxtbox.Text, mytxtbox.Text, Lxtxtbox.Text, Lytxtbox.Text);
                Thread RunThread = new Thread(RunModel);
                RunThread.Start(v);
                while (!RunThread.IsAlive) ;
                glControl1.Enabled = true;
            }
            else if (button1.Text == "Pause")
            {
                button1.Text = "Continue";
                pause = true;
            }
            else if (button1.Text == "Continue")
            {
                button1.Text = "Pause";
                pause = false;
            }
            
            
        }


        // Refresh Objects from Form2
        public void TimeTextRefresh(TextBox txtbox, string text)
        {
            if (txtbox.InvokeRequired)
                txtbox.Invoke((MethodInvoker)delegate()
                {
                    TimeTextRefresh(txtbox, text);
                });
            else
                txtbox.Text = text;
        }


        public void RunModel(object param)
        {
            RunParameter v = (RunParameter)param;
            Model M = new Model(v);
            double t;

            Fields.SetThermalField(Convert.ToDouble(TotextBox.Text), Convert.ToDouble(RtextBox.Text), Convert.ToDouble(GxtextBox.Text), Convert.ToDouble(GytextBox.Text));
            Output.CleanFile("test.dat", "Results");

            Model.Al myAlloy = new Model.Al();

            Kinectics.ThermalParameter tmp = Kinectics.Param(new string[] { "0.0", "8.26e-6","8.18e-5" });

            Mesh = BuildCAVisual(M.CA, v.mx, v.my, v.Lx, v.Ly);
            InitializePolygons();
            InitializedC();
            for (t = 0; t <= v.tmax + v.dtmax; t += v.dtmax)
            {
                
                tmp = Kinectics.Param( new string[] { Convert.ToString(Alloys.Tliq - Fields.ThermalField(t,0.0,0.0)), Convert.ToString(Alloys.a2), Convert.ToString(Alloys.a3) });
                M.Run(t);

                myAlloy.RefreshEsi(Fields.R, Fields.ThermalField(t,0.0,0.0), myAlloy.getA(M.CA, v), Convert.ToDouble(Kinectics.V(new object[] { tmp })), v.dtmax);

                if (Model.CheckTimeInterval(t, v.dtmax, Convert.ToDouble(rtitxtbox.Text)))
                {
                    TimeTextRefresh( textBox2, String.Format("{0:0.00}", t) );
                    TimeTextRefresh( textBox1, String.Format("{0:0.0000}", myAlloy.esi));

                    Output.WriteLineinFile(t, new double[] { myAlloy.esi }, "test.dat", "Results");

                    buildingVisual = true;
                    Visualization.RefreshStates(M.CA, state, v.mx, v.my, ref started);
                    Visualization.RefreshPolygons(M.CA, Mesh, Polygons, state, v.mx, v.my, v.Lx, v.Ly, dxGL);
                    Visualization.RefreshdC(M.CA, state, dC, v.mx, v.my, v.Lx, v.Ly, dxGL);
                    buildingVisual = false;
                }

                while (pause == true) ;
            }


            TimeTextRefresh( textBox2, M.Test());
        }


        // OpenGL Functions
        private void glControl1_Load(object sender, EventArgs e)
        {
            GLloaded = true;


            GL.ClearColor(Color.White);
            SetupViewPort();

            Application.Idle += Application_Idle;
        }

        void Application_Idle(object sendr, EventArgs e)
        {
            while (glControl1.IsIdle)
            {
                Render();
            }
        }

        private void Render()
        {
            if (!GLloaded)
                return;

            if (buildingVisual == false)
            {
                GL.Clear(ClearBufferMask.ColorBufferBit | ClearBufferMask.DepthBufferBit);

                GL.MatrixMode(MatrixMode.Modelview);
                GL.LoadIdentity();

                if (started)
                {
                    if (checkBox3.Checked == true)
                        Visualization.PointsVisualization(Mesh, mxp, myp);
                    if (checkBox4.Checked == true)
                        Visualization.SquaresVisualization(Mesh, mxp, myp, dxGL);

                        if (checkBox1.Checked == true)
                            Visualization.CAstateVisualization(Mesh, state, mxp, myp, dxGL);
                        if (checkBox2.Checked == true)
                            Visualization.CAPolygonVisualization(Polygons, state, mxp, myp);
                        if (checkBox5.Checked == true)
                            Visualization.dCVisualization(Mesh, dC, state, mxp, myp);

                
                    else
                    {

                    }
                }


                glControl1.SwapBuffers();
            }
            else
            {

            }
        }

        
        private void SetupViewPort()
        {
            int w = glControl1.Width;
            int h = glControl1.Height;
            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();
            GL.Ortho(0, w, 0, h, -1, 1);
            GL.Viewport(0, 0, w, h);
        }

        private void glControl1_Resize(object sender, EventArgs e)
        {
            if (!GLloaded)
                return;

        }

        private void glControl1_Paint(object sender, PaintEventArgs e)
        {
            Render();
        }

        private Vector2[,] BuildCAVisual(RappazGandinCACell[,] CA, int mx, int my, double Lx, double Ly)
        {
            int w = glControl1.Width;
            int h = glControl1.Height;
            Vector2[,] Grid = new Vector2[my, mx];
            mxp = mx;
            myp = my;
            Lyp = Ly;
            Lxp = Lx;
            for (int i = 0; i < my; i++)
            {
                for (int j = 0; j < mx; j++)
                {
                    Grid[i, j] = new Vector2();
                }
            }
            state = new int[my, mx];

            Vector2 Origin = new Vector2();
            double dxCA = Lx / mx;
            

            double ratioCA = Lx / Ly;
            double ratioGL = ((double) w) /((double) h);

            if (ratioGL > ratioCA)
            {
                dxGL = ((double) h) / my;
                Origin.v = new double[] { 0.5*(w - Lx/Ly*h) , 0.0 };
                LyVisual = h;
                LxVisual = Lx / Ly * h;
            }
            else
            {
                dxGL = ((double) w) / mx;
                Origin.v = new double[] { 0.0,  0.5 * (h - Ly / Lx * w) };
                LyVisual = w*Ly/Lx;
                LxVisual = h;
            }

            for (int i = 0; i < my; i++)
            {
                for (int j = 0; j < mx; j++)
                {
                    Grid[i,j].x = Origin.x + ((double) j + 0.5) * dxGL;
                    Grid[i,j].y = Origin.y + LyVisual - ((double) i + 0.5)*dxGL;
                }
            }

            return Grid;

        }

        void InitializePolygons()
        {
            Polygons = new Vector2[myp, mxp][];
            for( int i  = 0 ; i < myp ; i++){
                for (int j = 0; j < mxp; j++)
                {
                    Polygons[i,j] = new Vector2[4];

                    for (int k = 0; k < 4; k++)
                    {
                        Polygons[i, j][k] = new Vector2();
                    }

                }
            }
        }
        void InitializedC()
        {
            dC = new Vector2[myp, mxp];
            for (int i = 0; i < myp; i++)
            {
                for (int j = 0; j < mxp; j++)
                {
                    dC[i, j] = new Vector2();
                }
            }
        }

        // End of OpenGL Functions

        private void Form2_Load(object sender, EventArgs e)
        {
            Type mytype = typeof(Kinectics);
            MethodInfo[] NestedTypes = mytype.GetMethods();
            foreach (MethodInfo ty in NestedTypes)
            {
                if( ty.Name != "V" && ty.Name != "Set" && ty.Name != "GetHashCode" && ty.Name != "ToString" && ty.Name != "GetType" && ty.Name != "Equals")
                    comboBox1.Items.Add(ty.Name);
            }

            RtextBox.Text = Convert.ToString((Convert.ToDouble(TfinaltextBox.Text) - Convert.ToDouble(TotextBox.Text)) / Convert.ToDouble(tmaxtxtbox.Text));
            GxtextBox.Text = Convert.ToString((Convert.ToDouble(TrighttextBox.Text) - Convert.ToDouble(TlefttextBox.Text)) / Convert.ToDouble(Lxtxtbox.Text));
            GytextBox.Text = Convert.ToString((Convert.ToDouble(TuptextBox.Text) - Convert.ToDouble(TdowntextBox.Text)) / Convert.ToDouble(Lytxtbox.Text));

        }

        private void comboBox1_TextChanged(object sender, EventArgs e)
        {
            Kinectics.Set(comboBox1.Text);
        }

        private void button2_Click(object sender, EventArgs e)
        {

        }

        private void quitbutton_Click(object sender, EventArgs e)
        {
            this.Close();
        }

        private void TfinaltextBox_TextChanged(object sender, EventArgs e)
        {
            
        }

        private void TfinaltextBox_Leave(object sender, EventArgs e)
        {
            RtextBox.Text = Convert.ToString((Convert.ToDouble(TfinaltextBox.Text) - Convert.ToDouble(TotextBox.Text)) / Convert.ToDouble(tmaxtxtbox.Text));
            GxtextBox.Text = Convert.ToString((Convert.ToDouble(TrighttextBox.Text) - Convert.ToDouble(TlefttextBox.Text)) / Convert.ToDouble(Lxtxtbox.Text));
            GytextBox.Text = Convert.ToString((Convert.ToDouble(TuptextBox.Text) - Convert.ToDouble(TdowntextBox.Text)) / Convert.ToDouble(Lytxtbox.Text));
        }
    }
}
