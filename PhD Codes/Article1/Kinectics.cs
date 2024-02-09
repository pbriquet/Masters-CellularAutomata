using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text;
using System.Reflection;

namespace Article1
{
    static class Kinectics
    {

        static MethodInfo func;
        static bool setted = false;

        public static ThermalParameter T;
        // Funções disponíveis para o crescimento
        // Thermal
        public struct ThermalParameter
        {
            public double dT;
            public double a2;
            public double a3;
        };

        public static ThermalParameter Param( string[] par){
            ThermalParameter tmp;
            tmp.dT = Convert.ToDouble(par[0]); // First Parameters is dT
            tmp.a2 = Convert.ToDouble(par[1]); // Second Parameter is a2
            tmp.a3 = Convert.ToDouble(par[2]); // Third Parameter is a3

            T = tmp;
            return tmp;

        }
        public static double Thermal(ThermalParameter P){
            if (P.dT > 0.0)
                return P.a2 * P.dT * P.dT + P.a3 * P.dT * P.dT * P.dT;
            else
                return 0.0;
        }

        // Solutal
        public class SolutalParameters
        {
            public double Cd;
            public double Co;
            public double kpart;

            public void SolutalVariables(double iCd)
            {
                Cd = iCd;
            }
        };
        public static double Solutal(SolutalParameters P){
            return (P.Cd - P.Co);
        }


        // Configura qual é o tipo de cinética a ser utilizada
        public static void Set(string MethodName)
        {
            Type t = (typeof(Kinectics)); // Pega o nome da classe que contém o método
            setted = true;  // Indica que a cinética foi configurada.
            func = t.GetMethod(MethodName);     // func apontará para o função com o nome escrito no MethodName
             
        }

        // Método geral para cinética de crescimento
        public static object V(object[] methodArgs)
        {
            if (setted)
            {
                return func.Invoke(null, methodArgs);       // retorna a função escolhida em Set, com os parâmetros generalizadas methodArgs
            }
            else
            {
                throw new Exception("Kinectics not setted.");
            }
        }

        

    }
}
