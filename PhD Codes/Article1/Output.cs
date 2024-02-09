using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Reflection;

namespace Article1
{
    class Output
    {
        static string separator = ";";
        // Create a Directory in .exe Folder with string Name, and returns the Path.

        public static bool IsFileReady(String sFilename)
        {
            // If the file can be opened for exclusive access it means that the file
            // is no longer locked by another process.
            try
            {
                using (FileStream inputStream = File.Open(sFilename, FileMode.Open, FileAccess.Read, FileShare.None))
                {
                    if (inputStream.Length > 0)
                    {
                        return true;
                    }
                    else
                    {
                        return true;
                    }

                }
            }
            catch (Exception)
            {
                return false;
            }
        }

        public static string CreateOrGetFolder(string name)
        {
            // This is the full directory and exe name
            //String fullAppName = Assembly.GetExecutingAssembly().GetName().CodeBase;
    
            // This strips off the exe name
            //String fullAppPath = Path.GetDirectoryName(fullAppName);
            
            String fullAppPath = "C://Article//"; // Temporary Solution For Results
            // Combine the Path with new folder name
            string newPath = System.IO.Path.Combine(fullAppPath, name);
            // Correct Create Directory Argument
            Uri uri = new Uri(newPath);
            // Create new diretory with folder name
            if(!Directory.Exists(uri.LocalPath))
                Directory.CreateDirectory(uri.LocalPath);
            return newPath;

        }

        public static string CreateOrGetFile(string name, string folder){
            string Path = CreateOrGetFolder(folder);
            Path = System.IO.Path.Combine(Path, name);
            Uri uri = new Uri(Path);
            if (!System.IO.File.Exists(Path))
            {
                File.Create(uri.LocalPath);
            }

            return Path;
        }



        public static void WriteLineinFile(double x, double[] y, string filename, string foldername)
        {
            
            Uri uri = new Uri(CreateOrGetFile(filename,foldername));

                while(!IsFileReady(uri.LocalPath));

                using (System.IO.StreamWriter file = new System.IO.StreamWriter(uri.LocalPath, true))
                {
                    file.Write(Convert.ToString(x) + separator);
                    foreach (double yi in y){
                        file.Write(Convert.ToString(yi) + separator);
                    }
                    file.WriteLine("");
                    file.Close();

                }

        }

        public static void CleanFile(string filename, string folder)
        {
            Uri path = new Uri(CreateOrGetFile(filename, CreateOrGetFolder(folder)));

            while (!IsFileReady(path.LocalPath)) ;
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(path.LocalPath))
            {
                file.Write("");
                file.Close();
            }
        }

        public static void ChangeSeparator(string newsp)
        {
            separator = newsp;
        }

    }
}
