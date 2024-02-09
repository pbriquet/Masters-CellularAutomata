using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Threading;
using System.Reflection;
using System.Linq.Expressions;

namespace Article1
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        Form2 form2 = new Form2();

        public void startForm2()
        {
            Application.Run(form2);
        }

        private void button1_Click(object sender, EventArgs e)
        {
            
            Thread frm2 = new Thread(new ThreadStart( startForm2 ) );
            frm2.Start();

        }

        private void quitbutton_Click(object sender, EventArgs e)
        {
            Application.Exit();
        }


    }
}
