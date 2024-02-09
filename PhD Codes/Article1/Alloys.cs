using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Article1
{

    class Alloys
    {

        public static double a2 = 8.26e-6;
        public static double a3 = 8.18e-5;
        public static double Tliq = 891.0;
        /*
        public struct PropertiesParameters
        {
            double T;
            double C;
        }

        //delegate double Property(PropertiesParameters p);

        interface IAlloy
        {

            double Co();
            string Master();
            string Solute();

            double ps(PropertiesParameters p);
            double pl(PropertiesParameters p);
            
            
        }

        public class AlSi7 : IAlloy
        {

            public double Co(){ return 7.0; }
            public string Master() { return "Al"; }
            public string Solute() { return "Si"; }

            public double ps(PropertiesParameters p) { return 150.0; }
            public double pl(PropertiesParameters p) { return 160.0; }
            
            
        }

        public abstract class AAlloy
        {
            // General
            protected abstract double Co { public get; set; }
            protected abstract string Master { public get; set; }
            protected abstract string Solute { public get; set; }
            // Probable Parameters
            public double T { get; set; }
            public double C { get; set; }

            // General Properties
            public abstract double ps();
            public abstract double pl();

            // Thermodynamic Properties
            public abstract double Cp();
            public abstract double Lf();

            // Kinectic Properties

            // Thermal Properties


            // Solutal Properties
        }

        public static class AlSi7p : AAlloy
        {
            public AlSi7p() { this.Co = 7.0; this.Master = "Al"; this.Solute = "Si"; }

            public override double ps() { return 1500.0; } 
        }
         * 
         */
    }
  
}
