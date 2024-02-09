namespace Article1
{
    partial class Form2
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.button1 = new System.Windows.Forms.Button();
            this.glControl1 = new OpenTK.GLControl();
            this.label1 = new System.Windows.Forms.Label();
            this.textBox2 = new System.Windows.Forms.TextBox();
            this.label2 = new System.Windows.Forms.Label();
            this.tmaxtxtbox = new System.Windows.Forms.TextBox();
            this.dtmaxtxtbox = new System.Windows.Forms.TextBox();
            this.label3 = new System.Windows.Forms.Label();
            this.rtitxtbox = new System.Windows.Forms.TextBox();
            this.label4 = new System.Windows.Forms.Label();
            this.mxtxtbox = new System.Windows.Forms.TextBox();
            this.mytxtbox = new System.Windows.Forms.TextBox();
            this.Lxtxtbox = new System.Windows.Forms.TextBox();
            this.Lytxtbox = new System.Windows.Forms.TextBox();
            this.label5 = new System.Windows.Forms.Label();
            this.label6 = new System.Windows.Forms.Label();
            this.label7 = new System.Windows.Forms.Label();
            this.label8 = new System.Windows.Forms.Label();
            this.label9 = new System.Windows.Forms.Label();
            this.comboBox1 = new System.Windows.Forms.ComboBox();
            this.checkBox1 = new System.Windows.Forms.CheckBox();
            this.checkBox2 = new System.Windows.Forms.CheckBox();
            this.checkBox3 = new System.Windows.Forms.CheckBox();
            this.label10 = new System.Windows.Forms.Label();
            this.checkBox4 = new System.Windows.Forms.CheckBox();
            this.checkBox5 = new System.Windows.Forms.CheckBox();
            this.label11 = new System.Windows.Forms.Label();
            this.textBox1 = new System.Windows.Forms.TextBox();
            this.quitbutton = new System.Windows.Forms.Button();
            this.TotextBox = new System.Windows.Forms.TextBox();
            this.RtextBox = new System.Windows.Forms.TextBox();
            this.GxtextBox = new System.Windows.Forms.TextBox();
            this.GytextBox = new System.Windows.Forms.TextBox();
            this.label12 = new System.Windows.Forms.Label();
            this.label13 = new System.Windows.Forms.Label();
            this.label14 = new System.Windows.Forms.Label();
            this.label15 = new System.Windows.Forms.Label();
            this.TfinaltextBox = new System.Windows.Forms.TextBox();
            this.label16 = new System.Windows.Forms.Label();
            this.label17 = new System.Windows.Forms.Label();
            this.TlefttextBox = new System.Windows.Forms.TextBox();
            this.label18 = new System.Windows.Forms.Label();
            this.TuptextBox = new System.Windows.Forms.TextBox();
            this.label19 = new System.Windows.Forms.Label();
            this.TdowntextBox = new System.Windows.Forms.TextBox();
            this.label20 = new System.Windows.Forms.Label();
            this.TrighttextBox = new System.Windows.Forms.TextBox();
            this.SuspendLayout();
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(98, 266);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(147, 53);
            this.button1.TabIndex = 1;
            this.button1.Text = "Run";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // glControl1
            // 
            this.glControl1.BackColor = System.Drawing.Color.Black;
            this.glControl1.Location = new System.Drawing.Point(485, 12);
            this.glControl1.Name = "glControl1";
            this.glControl1.Size = new System.Drawing.Size(700, 700);
            this.glControl1.TabIndex = 3;
            this.glControl1.VSync = false;
            this.glControl1.Load += new System.EventHandler(this.glControl1_Load);
            this.glControl1.Paint += new System.Windows.Forms.PaintEventHandler(this.glControl1_Paint);
            this.glControl1.Resize += new System.EventHandler(this.glControl1_Resize);
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(307, 43);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(30, 13);
            this.label1.TabIndex = 2;
            this.label1.Text = "Time";
            // 
            // textBox2
            // 
            this.textBox2.Enabled = false;
            this.textBox2.Location = new System.Drawing.Point(355, 40);
            this.textBox2.Name = "textBox2";
            this.textBox2.Size = new System.Drawing.Size(97, 20);
            this.textBox2.TabIndex = 4;
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(95, 41);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(29, 13);
            this.label2.TabIndex = 5;
            this.label2.Text = "tmax";
            // 
            // tmaxtxtbox
            // 
            this.tmaxtxtbox.Location = new System.Drawing.Point(130, 38);
            this.tmaxtxtbox.Name = "tmaxtxtbox";
            this.tmaxtxtbox.Size = new System.Drawing.Size(93, 20);
            this.tmaxtxtbox.TabIndex = 6;
            this.tmaxtxtbox.Text = "100";
            this.tmaxtxtbox.TextChanged += new System.EventHandler(this.TfinaltextBox_TextChanged);
            this.tmaxtxtbox.Leave += new System.EventHandler(this.TfinaltextBox_Leave);
            // 
            // dtmaxtxtbox
            // 
            this.dtmaxtxtbox.Location = new System.Drawing.Point(130, 64);
            this.dtmaxtxtbox.Name = "dtmaxtxtbox";
            this.dtmaxtxtbox.Size = new System.Drawing.Size(93, 20);
            this.dtmaxtxtbox.TabIndex = 8;
            this.dtmaxtxtbox.Text = "0.00001";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(95, 67);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(35, 13);
            this.label3.TabIndex = 7;
            this.label3.Text = "dtmax";
            // 
            // rtitxtbox
            // 
            this.rtitxtbox.Cursor = System.Windows.Forms.Cursors.PanNorth;
            this.rtitxtbox.Location = new System.Drawing.Point(130, 90);
            this.rtitxtbox.Name = "rtitxtbox";
            this.rtitxtbox.Size = new System.Drawing.Size(93, 20);
            this.rtitxtbox.TabIndex = 10;
            this.rtitxtbox.Text = "0.01";
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(16, 93);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(108, 13);
            this.label4.TabIndex = 9;
            this.label4.Text = "Refresh Time Interval";
            // 
            // mxtxtbox
            // 
            this.mxtxtbox.Location = new System.Drawing.Point(132, 128);
            this.mxtxtbox.Name = "mxtxtbox";
            this.mxtxtbox.Size = new System.Drawing.Size(90, 20);
            this.mxtxtbox.TabIndex = 11;
            this.mxtxtbox.Text = "31";
            // 
            // mytxtbox
            // 
            this.mytxtbox.Location = new System.Drawing.Point(133, 154);
            this.mytxtbox.Name = "mytxtbox";
            this.mytxtbox.Size = new System.Drawing.Size(90, 20);
            this.mytxtbox.TabIndex = 12;
            this.mytxtbox.Text = "31";
            // 
            // Lxtxtbox
            // 
            this.Lxtxtbox.Location = new System.Drawing.Point(133, 180);
            this.Lxtxtbox.Name = "Lxtxtbox";
            this.Lxtxtbox.Size = new System.Drawing.Size(90, 20);
            this.Lxtxtbox.TabIndex = 13;
            this.Lxtxtbox.Text = "0.1";
            // 
            // Lytxtbox
            // 
            this.Lytxtbox.Location = new System.Drawing.Point(133, 206);
            this.Lytxtbox.Name = "Lytxtbox";
            this.Lytxtbox.Size = new System.Drawing.Size(90, 20);
            this.Lytxtbox.TabIndex = 14;
            this.Lytxtbox.Text = "0.1";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(74, 131);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(20, 13);
            this.label5.TabIndex = 15;
            this.label5.Text = "mx";
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(74, 157);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(20, 13);
            this.label6.TabIndex = 16;
            this.label6.Text = "my";
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(74, 187);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(18, 13);
            this.label7.TabIndex = 17;
            this.label7.Text = "Lx";
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(74, 213);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(18, 13);
            this.label8.TabIndex = 18;
            this.label8.Text = "Ly";
            // 
            // label9
            // 
            this.label9.AutoSize = true;
            this.label9.Location = new System.Drawing.Point(267, 131);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(50, 13);
            this.label9.TabIndex = 20;
            this.label9.Text = "Kinectics";
            // 
            // comboBox1
            // 
            this.comboBox1.FormattingEnabled = true;
            this.comboBox1.Location = new System.Drawing.Point(331, 127);
            this.comboBox1.Name = "comboBox1";
            this.comboBox1.Size = new System.Drawing.Size(105, 21);
            this.comboBox1.TabIndex = 21;
            this.comboBox1.TextChanged += new System.EventHandler(this.comboBox1_TextChanged);
            // 
            // checkBox1
            // 
            this.checkBox1.AutoSize = true;
            this.checkBox1.Checked = true;
            this.checkBox1.CheckState = System.Windows.Forms.CheckState.Checked;
            this.checkBox1.Location = new System.Drawing.Point(346, 190);
            this.checkBox1.Name = "checkBox1";
            this.checkBox1.Size = new System.Drawing.Size(68, 17);
            this.checkBox1.TabIndex = 22;
            this.checkBox1.Text = "CA State";
            this.checkBox1.UseVisualStyleBackColor = true;
            // 
            // checkBox2
            // 
            this.checkBox2.AutoSize = true;
            this.checkBox2.Checked = true;
            this.checkBox2.CheckState = System.Windows.Forms.CheckState.Checked;
            this.checkBox2.Location = new System.Drawing.Point(346, 213);
            this.checkBox2.Name = "checkBox2";
            this.checkBox2.Size = new System.Drawing.Size(86, 17);
            this.checkBox2.TabIndex = 23;
            this.checkBox2.Text = "CA Polygons";
            this.checkBox2.UseVisualStyleBackColor = true;
            // 
            // checkBox3
            // 
            this.checkBox3.AutoSize = true;
            this.checkBox3.Checked = true;
            this.checkBox3.CheckState = System.Windows.Forms.CheckState.Checked;
            this.checkBox3.Location = new System.Drawing.Point(353, 284);
            this.checkBox3.Name = "checkBox3";
            this.checkBox3.Size = new System.Drawing.Size(55, 17);
            this.checkBox3.TabIndex = 24;
            this.checkBox3.Text = "Points";
            this.checkBox3.UseVisualStyleBackColor = true;
            // 
            // label10
            // 
            this.label10.AutoSize = true;
            this.label10.Location = new System.Drawing.Point(343, 266);
            this.label10.Name = "label10";
            this.label10.Size = new System.Drawing.Size(26, 13);
            this.label10.TabIndex = 25;
            this.label10.Text = "Grid";
            // 
            // checkBox4
            // 
            this.checkBox4.AutoSize = true;
            this.checkBox4.Checked = true;
            this.checkBox4.CheckState = System.Windows.Forms.CheckState.Checked;
            this.checkBox4.Location = new System.Drawing.Point(352, 306);
            this.checkBox4.Name = "checkBox4";
            this.checkBox4.Size = new System.Drawing.Size(65, 17);
            this.checkBox4.TabIndex = 26;
            this.checkBox4.Text = "Squares";
            this.checkBox4.UseVisualStyleBackColor = true;
            // 
            // checkBox5
            // 
            this.checkBox5.AutoSize = true;
            this.checkBox5.Checked = true;
            this.checkBox5.CheckState = System.Windows.Forms.CheckState.Checked;
            this.checkBox5.Location = new System.Drawing.Point(346, 236);
            this.checkBox5.Name = "checkBox5";
            this.checkBox5.Size = new System.Drawing.Size(56, 17);
            this.checkBox5.TabIndex = 27;
            this.checkBox5.Text = "CA dC";
            this.checkBox5.UseVisualStyleBackColor = true;
            // 
            // label11
            // 
            this.label11.AutoSize = true;
            this.label11.Location = new System.Drawing.Point(303, 72);
            this.label11.Name = "label11";
            this.label11.Size = new System.Drawing.Size(20, 13);
            this.label11.TabIndex = 28;
            this.label11.Text = "esi";
            // 
            // textBox1
            // 
            this.textBox1.Enabled = false;
            this.textBox1.Location = new System.Drawing.Point(345, 70);
            this.textBox1.Name = "textBox1";
            this.textBox1.Size = new System.Drawing.Size(90, 20);
            this.textBox1.TabIndex = 29;
            // 
            // quitbutton
            // 
            this.quitbutton.Location = new System.Drawing.Point(95, 342);
            this.quitbutton.Name = "quitbutton";
            this.quitbutton.Size = new System.Drawing.Size(149, 50);
            this.quitbutton.TabIndex = 30;
            this.quitbutton.Text = "Quit";
            this.quitbutton.UseVisualStyleBackColor = true;
            this.quitbutton.Click += new System.EventHandler(this.quitbutton_Click);
            // 
            // TotextBox
            // 
            this.TotextBox.Location = new System.Drawing.Point(124, 450);
            this.TotextBox.Name = "TotextBox";
            this.TotextBox.Size = new System.Drawing.Size(59, 20);
            this.TotextBox.TabIndex = 31;
            this.TotextBox.Text = "881";
            this.TotextBox.Leave += new System.EventHandler(this.TfinaltextBox_Leave);
            // 
            // RtextBox
            // 
            this.RtextBox.Location = new System.Drawing.Point(351, 450);
            this.RtextBox.Name = "RtextBox";
            this.RtextBox.Size = new System.Drawing.Size(80, 20);
            this.RtextBox.TabIndex = 32;
            // 
            // GxtextBox
            // 
            this.GxtextBox.Location = new System.Drawing.Point(351, 476);
            this.GxtextBox.Name = "GxtextBox";
            this.GxtextBox.Size = new System.Drawing.Size(80, 20);
            this.GxtextBox.TabIndex = 33;
            // 
            // GytextBox
            // 
            this.GytextBox.Location = new System.Drawing.Point(352, 502);
            this.GytextBox.Name = "GytextBox";
            this.GytextBox.Size = new System.Drawing.Size(80, 20);
            this.GytextBox.TabIndex = 34;
            // 
            // label12
            // 
            this.label12.AutoSize = true;
            this.label12.Location = new System.Drawing.Point(98, 453);
            this.label12.Name = "label12";
            this.label12.Size = new System.Drawing.Size(20, 13);
            this.label12.TabIndex = 35;
            this.label12.Text = "To";
            // 
            // label13
            // 
            this.label13.AutoSize = true;
            this.label13.Location = new System.Drawing.Point(325, 453);
            this.label13.Name = "label13";
            this.label13.Size = new System.Drawing.Size(15, 13);
            this.label13.TabIndex = 36;
            this.label13.Text = "R";
            // 
            // label14
            // 
            this.label14.AutoSize = true;
            this.label14.Location = new System.Drawing.Point(325, 479);
            this.label14.Name = "label14";
            this.label14.Size = new System.Drawing.Size(20, 13);
            this.label14.TabIndex = 37;
            this.label14.Text = "Gx";
            // 
            // label15
            // 
            this.label15.AutoSize = true;
            this.label15.Location = new System.Drawing.Point(326, 505);
            this.label15.Name = "label15";
            this.label15.Size = new System.Drawing.Size(20, 13);
            this.label15.TabIndex = 38;
            this.label15.Text = "Gy";
            // 
            // TfinaltextBox
            // 
            this.TfinaltextBox.Location = new System.Drawing.Point(228, 450);
            this.TfinaltextBox.Name = "TfinaltextBox";
            this.TfinaltextBox.Size = new System.Drawing.Size(55, 20);
            this.TfinaltextBox.TabIndex = 39;
            this.TfinaltextBox.Text = "881";
            this.TfinaltextBox.Leave += new System.EventHandler(this.TfinaltextBox_Leave);
            // 
            // label16
            // 
            this.label16.AutoSize = true;
            this.label16.Location = new System.Drawing.Point(189, 453);
            this.label16.Name = "label16";
            this.label16.Size = new System.Drawing.Size(33, 13);
            this.label16.TabIndex = 40;
            this.label16.Text = "Tfinal";
            // 
            // label17
            // 
            this.label17.AutoSize = true;
            this.label17.Location = new System.Drawing.Point(85, 479);
            this.label17.Name = "label17";
            this.label17.Size = new System.Drawing.Size(28, 13);
            this.label17.TabIndex = 42;
            this.label17.Text = "Tleft";
            // 
            // TlefttextBox
            // 
            this.TlefttextBox.Location = new System.Drawing.Point(124, 476);
            this.TlefttextBox.Name = "TlefttextBox";
            this.TlefttextBox.Size = new System.Drawing.Size(58, 20);
            this.TlefttextBox.TabIndex = 41;
            this.TlefttextBox.Text = "891.0";
            this.TlefttextBox.Leave += new System.EventHandler(this.TfinaltextBox_Leave);
            // 
            // label18
            // 
            this.label18.AutoSize = true;
            this.label18.Location = new System.Drawing.Point(85, 505);
            this.label18.Name = "label18";
            this.label18.Size = new System.Drawing.Size(26, 13);
            this.label18.TabIndex = 44;
            this.label18.Text = "Tup";
            // 
            // TuptextBox
            // 
            this.TuptextBox.Location = new System.Drawing.Point(124, 502);
            this.TuptextBox.Name = "TuptextBox";
            this.TuptextBox.Size = new System.Drawing.Size(58, 20);
            this.TuptextBox.TabIndex = 43;
            this.TuptextBox.Text = "891.0";
            this.TuptextBox.Leave += new System.EventHandler(this.TfinaltextBox_Leave);
            // 
            // label19
            // 
            this.label19.AutoSize = true;
            this.label19.Location = new System.Drawing.Point(189, 505);
            this.label19.Name = "label19";
            this.label19.Size = new System.Drawing.Size(40, 13);
            this.label19.TabIndex = 48;
            this.label19.Text = "Tdown";
            // 
            // TdowntextBox
            // 
            this.TdowntextBox.Location = new System.Drawing.Point(228, 502);
            this.TdowntextBox.Name = "TdowntextBox";
            this.TdowntextBox.Size = new System.Drawing.Size(58, 20);
            this.TdowntextBox.TabIndex = 47;
            this.TdowntextBox.Text = "891.0";
            this.TdowntextBox.Leave += new System.EventHandler(this.TfinaltextBox_Leave);
            // 
            // label20
            // 
            this.label20.AutoSize = true;
            this.label20.Location = new System.Drawing.Point(189, 479);
            this.label20.Name = "label20";
            this.label20.Size = new System.Drawing.Size(34, 13);
            this.label20.TabIndex = 46;
            this.label20.Text = "Tright";
            // 
            // TrighttextBox
            // 
            this.TrighttextBox.Location = new System.Drawing.Point(228, 476);
            this.TrighttextBox.Name = "TrighttextBox";
            this.TrighttextBox.Size = new System.Drawing.Size(58, 20);
            this.TrighttextBox.TabIndex = 45;
            this.TrighttextBox.Text = "891.0";
            this.TrighttextBox.Leave += new System.EventHandler(this.TfinaltextBox_Leave);
            // 
            // Form2
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1399, 750);
            this.Controls.Add(this.label19);
            this.Controls.Add(this.TdowntextBox);
            this.Controls.Add(this.label20);
            this.Controls.Add(this.TrighttextBox);
            this.Controls.Add(this.label18);
            this.Controls.Add(this.TuptextBox);
            this.Controls.Add(this.label17);
            this.Controls.Add(this.TlefttextBox);
            this.Controls.Add(this.label16);
            this.Controls.Add(this.TfinaltextBox);
            this.Controls.Add(this.label15);
            this.Controls.Add(this.label14);
            this.Controls.Add(this.label13);
            this.Controls.Add(this.label12);
            this.Controls.Add(this.GytextBox);
            this.Controls.Add(this.GxtextBox);
            this.Controls.Add(this.RtextBox);
            this.Controls.Add(this.TotextBox);
            this.Controls.Add(this.quitbutton);
            this.Controls.Add(this.textBox1);
            this.Controls.Add(this.label11);
            this.Controls.Add(this.checkBox5);
            this.Controls.Add(this.checkBox4);
            this.Controls.Add(this.label10);
            this.Controls.Add(this.checkBox3);
            this.Controls.Add(this.checkBox2);
            this.Controls.Add(this.checkBox1);
            this.Controls.Add(this.comboBox1);
            this.Controls.Add(this.label9);
            this.Controls.Add(this.label8);
            this.Controls.Add(this.label7);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.Lytxtbox);
            this.Controls.Add(this.Lxtxtbox);
            this.Controls.Add(this.mytxtbox);
            this.Controls.Add(this.mxtxtbox);
            this.Controls.Add(this.rtitxtbox);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.dtmaxtxtbox);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.tmaxtxtbox);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.textBox2);
            this.Controls.Add(this.glControl1);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.button1);
            this.Name = "Form2";
            this.Text = "Form2";
            this.Load += new System.EventHandler(this.Form2_Load);
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button button1;
        private OpenTK.GLControl glControl1;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.TextBox textBox2;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.TextBox tmaxtxtbox;
        private System.Windows.Forms.TextBox dtmaxtxtbox;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.TextBox rtitxtbox;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.TextBox mxtxtbox;
        private System.Windows.Forms.TextBox mytxtbox;
        private System.Windows.Forms.TextBox Lxtxtbox;
        private System.Windows.Forms.TextBox Lytxtbox;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.Label label9;
        private System.Windows.Forms.ComboBox comboBox1;
        private System.Windows.Forms.CheckBox checkBox1;
        private System.Windows.Forms.CheckBox checkBox2;
        private System.Windows.Forms.CheckBox checkBox3;
        private System.Windows.Forms.Label label10;
        private System.Windows.Forms.CheckBox checkBox4;
        private System.Windows.Forms.CheckBox checkBox5;
        private System.Windows.Forms.Label label11;
        private System.Windows.Forms.TextBox textBox1;
        private System.Windows.Forms.Button quitbutton;
        private System.Windows.Forms.TextBox TotextBox;
        private System.Windows.Forms.TextBox RtextBox;
        private System.Windows.Forms.TextBox GxtextBox;
        private System.Windows.Forms.TextBox GytextBox;
        private System.Windows.Forms.Label label12;
        private System.Windows.Forms.Label label13;
        private System.Windows.Forms.Label label14;
        private System.Windows.Forms.Label label15;
        private System.Windows.Forms.TextBox TfinaltextBox;
        private System.Windows.Forms.Label label16;
        private System.Windows.Forms.Label label17;
        private System.Windows.Forms.TextBox TlefttextBox;
        private System.Windows.Forms.Label label18;
        private System.Windows.Forms.TextBox TuptextBox;
        private System.Windows.Forms.Label label19;
        private System.Windows.Forms.TextBox TdowntextBox;
        private System.Windows.Forms.Label label20;
        private System.Windows.Forms.TextBox TrighttextBox;

    }
}