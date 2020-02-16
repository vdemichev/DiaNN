using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;

using EasyTabs;

namespace GUI
{
    static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        public static string pipeline_to_open = ""; 

        [STAThread]
        static void Main()
        {
            var culture = System.Globalization.CultureInfo.CreateSpecificCulture("en-US");
            System.Globalization.CultureInfo.DefaultThreadCurrentCulture = culture;
            System.Globalization.CultureInfo.DefaultThreadCurrentUICulture = culture;
            System.Threading.Thread.CurrentThread.CurrentCulture = culture;
            System.Threading.Thread.CurrentThread.CurrentUICulture = culture;

            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);

            AppContainer container = new AppContainer();

            string[] args = Environment.GetCommandLineArgs();
            if (args.Length > 1) if (args[1].Contains(".pipeline")) pipeline_to_open = args[1];

            // Add the initial Tab
            container.Tabs.Add(
                // Our First Tab created by default in the Application will have as content the Form1
                new TitleBarTab(container)
                {
                    Content = new Form1 { }
                }
            );
            pipeline_to_open = "";

            // Set initial tab the first one
            container.SelectedTabIndex = 0;

            // Create tabs and start application
            TitleBarTabsApplicationContext applicationContext = new TitleBarTabsApplicationContext();
            applicationContext.Start(container);

            Application.Run(applicationContext);
        }
    }
}
