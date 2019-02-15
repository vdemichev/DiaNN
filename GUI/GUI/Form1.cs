using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Windows.Forms;
using System.Linq;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;

namespace GUI
{
    [Serializable]
    public struct Settings
    {
        public string name;
        public string files_s, lib_s, out_s, out_gene_s, fasta_s, diann_s, learn_lib_s, out_lib_s, add_s;
        public decimal threads_i, log_i, missed_i, varmod_i, scan_i, pep_min_i, pep_max_i,
            pr_min_i, pr_max_i, fr_min_i, fr_max_i;
        public int grouping_i, protease_i;
        public bool batches_b, nn_b, unrelated_b, rtprof_b, protinf_b, iso_b, ifr_b,
            use_quant_b, use_lib_free_b, met_exc_b, carbamet_b, oxid_b, opt_training_b;
        public decimal prec_fdr_d, prot_fdr_d, mass_acc_d, mass_acc_ms1_d;
    }

    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            int cores = 0;
            foreach (var pr in new System.Management.ManagementObjectSearcher("Select * from Win32_Processor").Get())
                cores += int.Parse(pr["NumberOfCores"].ToString());
            ThreadsUpDown.Value = cores;
            ThreadsUpDown.Maximum = Environment.ProcessorCount;
            EnzymeCombo.SelectedIndex = 0;
            PGBox.SelectedIndex = 3;
        }

        bool running = false, finished = false, in_pipeline = false, finished_pipeline = false;
        Process process = new Process();
        List<Settings> Pipeline = new List<Settings>();
        List<string> PipNames = new List<string>();
        int pip_last_index = 0, processed = 0;
        
        private void SaveSettings(ref Settings S)
        {
            S.diann_s = ExeText.Text;
            S.files_s = FilesList.Text;
            S.lib_s = LibText.Text;
            S.threads_i = ThreadsUpDown.Value;
            S.log_i = VerboseUpDown.Value;
            S.out_s = OutputText.Text;
            S.out_gene_s = GeneRepBox.Text;
            S.prec_fdr_d = FDRUpDown.Value;
            S.prot_fdr_d = ProteinFDRUpDown.Value;
            S.fasta_s = FastaText.Text;
            S.use_lib_free_b = LibraryFreeBox.Checked;
            S.out_lib_s = OutputLibText.Text;
            S.learn_lib_s = LearnLibText.Text;
            S.protease_i = EnzymeCombo.SelectedIndex;
            S.missed_i = MissedCleavageUpDown.Value;
            S.pep_min_i = PepLenMin.Value;
            S.pep_max_i = PepLenMax.Value;
            S.pr_min_i = PrMzMin.Value;
            S.pr_max_i = PrMzMax.Value;
            S.fr_min_i = FrMzMin.Value;
            S.fr_max_i = FrMzMax.Value;
            S.met_exc_b = MethionineBox.Checked;
            S.carbamet_b = CysteineBox.Checked;
            S.varmod_i = VarModsUpDown.Value;
            S.oxid_b = OxidationBox.Checked;
            S.opt_training_b = OptimiseTrainingBox.Checked;
            S.scan_i = WindowUpDown.Value;
            S.mass_acc_d = MassAccUpDown.Value;
            S.mass_acc_ms1_d = MassAccMs1UpDown.Value;
            S.use_quant_b = UseQuantCheck.Checked;
            S.batches_b = BatchModeCheck.Checked;
            S.nn_b = NNCheck.Checked;
            S.unrelated_b = IndividualRunsCheck.Checked;
            S.rtprof_b = RTProfilingBox.Checked;
            S.protinf_b = ProtInfBox.Checked;
            S.iso_b = IsotopeBox.Checked;
            S.ifr_b = InterferenceBox.Checked;
            S.grouping_i = PGBox.SelectedIndex;
            S.add_s = OptionsText.Text;
        }

        private void LoadSettings(Settings S)
        {
            ExeText.Text = S.diann_s;
            FilesList.Text = S.files_s;
            LibText.Text = S.lib_s;
            ThreadsUpDown.Value = S.threads_i;
            VerboseUpDown.Value = S.log_i;
            OutputText.Text = S.out_s;
            GeneRepBox.Text = S.out_gene_s;
            FDRUpDown.Value = S.prec_fdr_d;
            ProteinFDRUpDown.Value = S.prot_fdr_d;
            FastaText.Text = S.fasta_s;
            LibraryFreeBox.Checked = S.use_lib_free_b;
            OutputLibText.Text = S.out_lib_s;
            LearnLibText.Text = S.learn_lib_s;
            EnzymeCombo.SelectedIndex = S.protease_i;
            MissedCleavageUpDown.Value = S.missed_i;
            PepLenMin.Value = S.pep_min_i;
            PepLenMax.Value = S.pep_max_i;
            PrMzMin.Value = S.pr_min_i;
            PrMzMax.Value = S.pr_max_i;
            FrMzMin.Value = S.fr_min_i;
            FrMzMax.Value = S.fr_max_i;
            MethionineBox.Checked = S.met_exc_b;
            CysteineBox.Checked = S.carbamet_b;
            VarModsUpDown.Value = S.varmod_i;
            OxidationBox.Checked = S.oxid_b;
            OptimiseTrainingBox.Checked = S.opt_training_b;
            WindowUpDown.Value = S.scan_i;
            MassAccUpDown.Value = S.mass_acc_d;
            MassAccMs1UpDown.Value = S.mass_acc_ms1_d;
            UseQuantCheck.Checked = S.use_quant_b;
            BatchModeCheck.Checked = S.batches_b;
            NNCheck.Checked = S.nn_b;
            IndividualRunsCheck.Checked = S.unrelated_b;
            RTProfilingBox.Checked = S.rtprof_b;
            ProtInfBox.Checked = S.protinf_b;
            IsotopeBox.Checked = S.iso_b;
            InterferenceBox.Checked = S.ifr_b;
            PGBox.SelectedIndex = S.grouping_i;
            OptionsText.Text = S.add_s;
        }

        private void RawDataButton_Click(object sender, EventArgs e)
        {
            OpenFileDialog rawDataDialog = new OpenFileDialog();
            rawDataDialog.Multiselect = true;
            rawDataDialog.Filter = "MS data files (*.raw, *.wiff, *.mzML, *.dia)|*.raw;*.wiff;*.mzML;*.dia|All files (*.*)|*.*";
            rawDataDialog.FilterIndex = 0;
            if (rawDataDialog.ShowDialog() == DialogResult.OK)
            {
                foreach (var str in rawDataDialog.FileNames)
                {
                    FilesList.AppendText(str + '\n');
                }
            }
        }

        private void LibButton_Click(object sender, EventArgs e)
        {
            OpenFileDialog libDialog = new OpenFileDialog();
            libDialog.Filter = "Text files (*.txt, *csv, *.tsv, *xls)|*.txt;*.csv;*.tsv;*.xls|All files (*.*)|*.*";
            libDialog.FilterIndex = 0;
            if (libDialog.ShowDialog() == DialogResult.OK)
                LibText.Text = libDialog.FileName;
        }

        private void OutputHandler(object sender, DataReceivedEventArgs outLine)
        {
            string line = (outLine?.Data?.ToString());
            if (line == "Finished") {
                finished = true;
                if (!running && (in_pipeline || finished_pipeline))
                    Invoke(new Action(() => PipelineList.Items[processed].BackColor = System.Drawing.SystemColors.Control));
            }
            Invoke(new Action(() => LogText.AppendText(line + Environment.NewLine)));
        }

        private void ExitedHandler(object sender, System.EventArgs e)
        {
            running = false;
            Invoke(new Action(() => LogText.Text += "DIA-NN exited" + Environment.NewLine));
            Invoke(new Action(() => StatusLabel.Text = "Finished"));
            Invoke(new Action(() => StatusIndicator.BackColor = System.Drawing.SystemColors.Control));
            if (in_pipeline)
            {
                bool failed = false;
                if (!finished)
                {
                    if (!process.HasExited) failed = true;
                    else if (process.ExitCode != 0) failed = true;
                    if (failed) Invoke(new Action(() => PipelineList.Items[processed].BackColor = System.Drawing.Color.LightSalmon));
                }
                if (!failed) Invoke(new Action(() => PipelineList.Items[processed].BackColor = System.Drawing.SystemColors.Control));
                processed++;
                if (processed < Pipeline.Count)
                {
                    process.Dispose();
                    Invoke(new Action(() => RunProcess(false, Pipeline[processed])));
                }
                else
                {
                    in_pipeline = false;
                    finished_pipeline = true;
                }
            }
        }

        private void RunProcess(bool convert, Settings S)
        {
            if (running) return;
            finished = finished_pipeline = false;

            if (in_pipeline) PipelineList.Items[processed].BackColor = System.Drawing.Color.GreenYellow;

            bool external = false;
            if (S.add_s.Length >= 1) if (S.add_s[0] == '>') external = true; // run some external tool; only "Additional options" will be past to the tool
            string opts;
            if (external) opts = S.add_s.Substring(1);
            else opts = S.add_s;

            if (!external && S.add_s.Contains("--convert")) convert = true;

            process = new Process();
            process.StartInfo.FileName = S.diann_s;
            process.StartInfo.CreateNoWindow = true;
            process.StartInfo.UseShellExecute = false;
            process.StartInfo.RedirectStandardOutput = true;
            process.StartInfo.RedirectStandardInput = true;
            process.EnableRaisingEvents = true;

            if (!external)
            {
                List<string> rawDataFiles = new List<string>(S.files_s.Split(new[] { '\n' }));
                rawDataFiles = rawDataFiles.Distinct(StringComparer.InvariantCultureIgnoreCase).ToList();
                foreach (var file in rawDataFiles)
                    if (file.Length >= 2) process.StartInfo.Arguments += " --f \"" + file + "\"";
            }
            if (convert)
            {
                process.StartInfo.Arguments += " --convert";
                process.StartInfo.Arguments += " --threads " + S.threads_i.ToString();
                process.StartInfo.Arguments += " --verbose " + S.log_i.ToString();
            }
            else
            {
                if (!external)
                {
                    process.StartInfo.Arguments += " --lib \"" + S.lib_s + "\"";
                    process.StartInfo.Arguments += " --threads " + S.threads_i.ToString();
                    process.StartInfo.Arguments += " --verbose " + S.log_i.ToString();
                    process.StartInfo.Arguments += " --out \"" + S.out_s + "\"";
                    process.StartInfo.Arguments += " --out-gene \"" + S.out_gene_s + "\"";
                    process.StartInfo.Arguments += " --qvalue " + Convert.ToString(0.01 * (double)S.prec_fdr_d);
                    process.StartInfo.Arguments += " --protein-qvalue " + Convert.ToString(0.01 * (double)S.prot_fdr_d);

                    if (S.fasta_s != "")
                    {
                        List<string> fastaFiles = new List<string>(S.fasta_s.Split(new[] { '\n' }));
                        foreach (var file in fastaFiles)
                            if (file.Length >= 2) process.StartInfo.Arguments += " --fasta \"" + file + "\"";
                        if (S.use_lib_free_b)
                        {
                            process.StartInfo.Arguments += " --fasta-search --out-lib \"" + S.out_lib_s + "\"";
                            if (S.learn_lib_s != "") process.StartInfo.Arguments += " --learn-lib \"" + S.learn_lib_s + "\"";

                            if (S.protease_i == 0) process.StartInfo.Arguments += " --cut-after KR";
                            if (S.protease_i == 1) process.StartInfo.Arguments += " --cut-after KR --no-cut-before P";
                            if (S.protease_i == 2) process.StartInfo.Arguments += " --cut-after K";
                            if (S.protease_i == 3) process.StartInfo.Arguments += " --cut-after K --no-cut-before P";
                            if (S.protease_i == 4) process.StartInfo.Arguments += " --cut-after FYWL --no-cut-before P";

                            process.StartInfo.Arguments += " --missed-cleavages " + S.missed_i.ToString();

                            process.StartInfo.Arguments += " --min-pep-len " + S.pep_min_i.ToString();
                            process.StartInfo.Arguments += " --max-pep-len " + S.pep_max_i.ToString();
                            process.StartInfo.Arguments += " --min-pr-mz " + S.pr_min_i.ToString();
                            process.StartInfo.Arguments += " --max-pr-mz " + S.pr_max_i.ToString();
                            process.StartInfo.Arguments += " --min-fr-mz " + S.fr_min_i.ToString();
                            process.StartInfo.Arguments += " --max-fr-mz " + S.fr_max_i.ToString();

                            if (S.opt_training_b) process.StartInfo.Arguments += " --min-fr-corr 0.9 --min-gen-fr 2";
                            if (S.met_exc_b) process.StartInfo.Arguments += " --met-excision";
                            if (S.carbamet_b) process.StartInfo.Arguments += " --unimod4";
                            if (S.varmod_i >= 1)
                            {
                                process.StartInfo.Arguments += " --var-mods " + S.varmod_i.ToString();
                                if (S.oxid_b) process.StartInfo.Arguments += " --unimod35";
                            }
                        }
                    }
                    else if (S.use_lib_free_b) process.StartInfo.Arguments += " --gen-spec-lib";

                    if ((int)S.scan_i > 0) process.StartInfo.Arguments += " --window " + S.scan_i.ToString();
                    if ((double)S.mass_acc_d > 0.0) process.StartInfo.Arguments += " --mass-acc " + S.mass_acc_d.ToString();
                    if ((double)S.mass_acc_ms1_d > 0.0) process.StartInfo.Arguments += " --mass-acc-ms1 " + S.mass_acc_ms1_d.ToString();

                    if (S.use_quant_b) process.StartInfo.Arguments += " --use-quant";
                    if (!S.batches_b) process.StartInfo.Arguments += " --no-batch-mode";
                    if (!S.nn_b) process.StartInfo.Arguments += " --no-nn";
                    if (S.unrelated_b) process.StartInfo.Arguments += " --individual-mass-acc --individual-windows";
                    if (S.rtprof_b) process.StartInfo.Arguments += " --rt-profiling";
                    if (!S.protinf_b) process.StartInfo.Arguments += " --no-prot-inf";
                    if (!S.iso_b) process.StartInfo.Arguments += " --no-isotopes";
                    if (!S.ifr_b) process.StartInfo.Arguments += " --int-removal 0";

                    if (S.grouping_i <= 2)
                    {
                        process.StartInfo.Arguments += " --pg-level " + S.grouping_i.ToString();
                        if (S.grouping_i == 2) process.StartInfo.Arguments += " --species-genes";
                    }
                }

                process.StartInfo.Arguments += " " + S.add_s;
            }

            process.OutputDataReceived += new DataReceivedEventHandler(OutputHandler);
            process.Exited += new EventHandler(ExitedHandler);

            if (process.StartInfo.Arguments.Length >= 32000 && !external)
            {
                LogText.Text += "Large number of files cannot be passed to the command line tool as command line arguments. " +
                    "A config file will therefore be created and referenced with a --cfg command. " +
                    "Alternatively, one can use --dir command to include all files in a directory." + Environment.NewLine;
                String fname = Directory.GetCurrentDirectory() + "\\DIA-NN-launch-cfg.txt";
                try
                {
                    File.WriteAllText(fname, process.StartInfo.Arguments);
                    process.StartInfo.Arguments = " --cfg \"" + fname + "\"";
                }
                catch (Exception ex) { LogText.Text += "ERROR: cannot create a config file, check if the folder with DIA-NN GUI is write-protected." + Environment.NewLine; }
            }

            try
            {
                running = process.Start();
            }
            catch (Exception ex) { running = false; }
            if (!running)
            {
                LogText.Text += "Failed to start " + S.diann_s + Environment.NewLine;
                if (in_pipeline)
                {
                    PipelineList.Items[processed].BackColor = System.Drawing.Color.LightSalmon;
                    processed++;
                    if (processed >= PipelineList.Items.Count) in_pipeline = false;
                    else
                    {
                        process.Dispose();
                        RunProcess(false, Pipeline[processed]);
                    }
                }
                return;
            }
            StatusLabel.Text = "Running...";
            StatusIndicator.BackColor = System.Drawing.Color.GreenYellow;
            LogText.Text += process.StartInfo.FileName + process.StartInfo.Arguments + Environment.NewLine;
            process.BeginOutputReadLine();
        }

        private void RunButton_Click(object sender, EventArgs e)
        {
            Settings S = new Settings();
            SaveSettings(ref S);
            process.Dispose();
            RunProcess(false, S);
        }

        private void ClearButton_Click(object sender, EventArgs e)
        {
            FilesList.Text = "";
        }

        private void StopButton_Click(object sender, EventArgs e)
        {
            if (!running) return;

            DialogResult dr;
            if (!in_pipeline) dr = MessageBox.Show("Stop DIA-NN? All progress will be lost.", "Stop", MessageBoxButtons.YesNo);
            else dr = MessageBox.Show("Cancel the current pipeline step? All progress will be lost.", "Stop", MessageBoxButtons.YesNo);
            switch (dr)
            {
                case DialogResult.Yes:
                    try
                    {
                        process.Kill();
                    }
                    catch (Exception ex) {}
                    return;
                case DialogResult.No:
                    return;
            }
        }

        private void LogClearButton_Click(object sender, EventArgs e)
        {
            LogText.Text = "";
        }

        private void ExeButton_Click(object sender, EventArgs e)
        {
            OpenFileDialog exeDialog = new OpenFileDialog();
            exeDialog.Filter = "exe files (*.exe)|*.exe|All files (*.*)|*.*";
            exeDialog.FilterIndex = 0;
            if (exeDialog.ShowDialog() == DialogResult.OK)
                ExeText.Text = exeDialog.FileName;
        }

        private void ConvertButton_Click(object sender, EventArgs e)
        {
            Settings S = new Settings();
            SaveSettings(ref S);
            process.Dispose();
            RunProcess(true, S);
        }

        private void OutputButton_Click(object sender, EventArgs e)
        {
            SaveFileDialog outDialog = new SaveFileDialog();
            outDialog.Filter = "Tab-separated files (*.tsv)|*.tsv|All files (*.*)|*.*";
            outDialog.FilterIndex = 0;
            outDialog.RestoreDirectory = true;
            if (outDialog.ShowDialog() == DialogResult.OK)
                OutputText.Text = outDialog.FileName;
        }

        private void FastaButton_Click(object sender, EventArgs e)
        {
            OpenFileDialog fastaDialog = new OpenFileDialog();
            fastaDialog.Multiselect = true;
            fastaDialog.Filter = "FASTA files (*.fasta, *fa)|*.fasta;*.fa|All files (*.*)|*.*";
            fastaDialog.FilterIndex = 0;
            if (fastaDialog.ShowDialog() == DialogResult.OK)
            {
                foreach (var str in fastaDialog.FileNames)
                {
                    FastaText.AppendText(str + '\n');
                }
            }
        }

        private void LearnLibButton_Click(object sender, EventArgs e)
        {
            OpenFileDialog libDialog = new OpenFileDialog();
            libDialog.Filter = "Text files (*.txt, *csv, *.tsv, *xls)|*.txt;*.csv;*.tsv;*.xls|All files (*.*)|*.*";
            libDialog.FilterIndex = 0;
            if (libDialog.ShowDialog() == DialogResult.OK)
                LearnLibText.Text = libDialog.FileName;
        }

        private void OutputLibButton_Click(object sender, EventArgs e)
        {
            SaveFileDialog outDialog = new SaveFileDialog();
            outDialog.Filter = "Tab-separated files (*.tsv)|*.tsv|All files (*.*)|*.*";
            outDialog.FilterIndex = 0;
            outDialog.RestoreDirectory = true;
            if (outDialog.ShowDialog() == DialogResult.OK)
                OutputLibText.Text = outDialog.FileName;
        }

        private void MassAccMs1UpDown_ValueChanged(object sender, EventArgs e)
        {
            if ((double)MassAccMs1UpDown.Value > 0.0 && !((double)MassAccUpDown.Value > 0.0)) MassAccUpDown.Value = (decimal)20.0;
        }

        private void MassAccUpDown_ValueChanged(object sender, EventArgs e)
        {
            if ((double)MassAccUpDown.Value > 0.0 && !((double)MassAccMs1UpDown.Value > 0.0)) MassAccMs1UpDown.Value = (decimal)20.0;
        }

        private void SaveLogButton_Click(object sender, EventArgs e)
        {
            SaveFileDialog outDialog = new SaveFileDialog();
            outDialog.Filter = "Text files (*.txt)|*.txt|All files (*.*)|*.*";
            outDialog.FilterIndex = 0;
            outDialog.RestoreDirectory = true;
            if (outDialog.ShowDialog() == DialogResult.OK)
                File.WriteAllText(outDialog.FileName, LogText.Text);
        }

        private void PepLenMin_ValueChanged(object sender, EventArgs e)
        {
            if (PepLenMin.Value > PepLenMax.Value) PepLenMax.Value = PepLenMin.Value;
        }

        private void PepLenMax_ValueChanged(object sender, EventArgs e)
        {
            if (PepLenMin.Value > PepLenMax.Value) PepLenMin.Value = PepLenMax.Value;
        }

        private void PrMzMin_ValueChanged(object sender, EventArgs e)
        {
            if (PrMzMin.Value > PrMzMax.Value) PrMzMax.Value = PrMzMin.Value;
        }

        private void PrMzMax_ValueChanged(object sender, EventArgs e)
        {
            if (PrMzMin.Value > PrMzMax.Value) PrMzMin.Value = PrMzMax.Value;
        }

        private void FrMzMin_ValueChanged(object sender, EventArgs e)
        {
            if (FrMzMin.Value > FrMzMax.Value) FrMzMax.Value = FrMzMin.Value;
        }

        private void FrMzMax_ValueChanged(object sender, EventArgs e)
        {
            if (FrMzMin.Value > FrMzMax.Value) FrMzMin.Value = FrMzMax.Value;
        }

        private void Form1_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (MessageBox.Show("Close DIA-NN? All progress will be lost.", "Close",
               MessageBoxButtons.YesNo, MessageBoxIcon.Question) == DialogResult.No)
            {
                e.Cancel = true;
            }
            else
            {
                try
                {
                    process.Kill();
                }
                catch (Exception ex) { }
            }
        }

        private void ClearFastaButton_Click(object sender, EventArgs e)
        {
            FastaText.Text = "";
        }

        private void OxidationBox_CheckedChanged(object sender, EventArgs e)
        {
            if (OxidationBox.Checked && VarModsUpDown.Value < 1)
                VarModsUpDown.Value = 1;
        }

        private void GeneRepButton_Click(object sender, EventArgs e)
        {
            SaveFileDialog outDialog = new SaveFileDialog();
            outDialog.Filter = "Tab-separated files (*.tsv)|*.tsv|All files (*.*)|*.*";
            outDialog.FilterIndex = 0;
            outDialog.RestoreDirectory = true;
            if (outDialog.ShowDialog() == DialogResult.OK)
                GeneRepBox.Text = outDialog.FileName;
        }

        private void PipAdd_Click(object sender, EventArgs e)
        {
            Settings S = new Settings();
            SaveSettings(ref S);
            if (PipelineList.SelectedIndices.Count == 0)
            {
                Pipeline.Add(S);
                PipNames.Add(PipName.Text);
                PipelineList.Items.Add(PipName.Text);
                PipelineList.Items[PipelineList.Items.Count - 1].Selected = true;
                PipelineList.Select();
            }
            else if (PipelineList.SelectedIndices.Count == 1)
                if (!in_pipeline || PipelineList.SelectedIndices[0] >= processed)
            {
                int i = PipelineList.SelectedIndices[0] + 1;
                Pipeline.Insert(i, S);
                PipNames.Insert(i, PipName.Text);
                PipelineList.Items.Insert(i, PipName.Text);
                PipelineList.Items[i - 1].Selected = false;
                PipelineList.Items[i].Selected = true;
                PipelineList.Select();
            }
            pip_last_index++;
            PipName.Text = "Step " + (pip_last_index + 1).ToString();
        }

        private void PipClear_Click(object sender, EventArgs e)
        {
            if (!in_pipeline && Pipeline.Count >= 1)
            {
                DialogResult dr = MessageBox.Show("Clear the entire pipeline?", "Clear", MessageBoxButtons.YesNo);
                switch (dr)
                {
                    case DialogResult.Yes:
                        try
                        {
                            Pipeline.Clear();
                            PipNames.Clear();
                            PipelineList.Items.Clear();
                            pip_last_index = 0;
                            PipName.Text = "Step 1";
                        }
                        catch (Exception ex) { }
                        return;
                    case DialogResult.No:
                        return;
                }
            }
        }

        private void PipRemove_Click(object sender, EventArgs e)
        {
            if (PipelineList.SelectedIndices.Count == 1)
            {
                int i = PipelineList.SelectedIndices[0];
                if (!in_pipeline || i > processed)
                {
                    PipelineList.Items.RemoveAt(i);
                    Pipeline.RemoveAt(i);
                    PipNames.RemoveAt(i);
                    pip_last_index--;
                }
            }
        }

        private void PipAbort_Click(object sender, EventArgs e)
        {
            if (!running || !in_pipeline) return;

            DialogResult dr = MessageBox.Show("Abort the pipeline? All progress (current step) will be lost.", "Abort", MessageBoxButtons.YesNo);
            switch (dr)
            {
                case DialogResult.Yes:
                    try
                    {
                        Invoke(new Action(() => PipelineList.Items[processed].BackColor = System.Drawing.Color.LightSalmon));
                        in_pipeline = false;
                        processed = 0;
                        process.Kill();
                    }
                    catch (Exception ex) { }
                    return;
                case DialogResult.No:
                    return;
            }
        }

        private void SavePipelineButton_Click(object sender, EventArgs e)
        {
            if (Pipeline.Count == 0) return;
            SaveFileDialog outDialog = new SaveFileDialog();
            outDialog.Filter = "Pipeline files (*.pipeline)|*.pipeline|All files (*.*)|*.*";
            outDialog.FilterIndex = 0;
            outDialog.RestoreDirectory = true;
            if (outDialog.ShowDialog() == DialogResult.OK)
            {

                int i;
                FileStream stream = File.Open(outDialog.FileName, FileMode.Create);
                var formatter = new BinaryFormatter();
                for (i = 0; i < Pipeline.Count; i++)
                {
                    var p = Pipeline[i];
                    p.name = PipNames[i];
                    formatter.Serialize(stream, p);
                }
                stream.Close();
            }
        }

        private void OpenPipelineButton_Click(object sender, EventArgs e)
        {
            OpenFileDialog openDialog = new OpenFileDialog();
            openDialog.Multiselect = false;
            openDialog.Filter = "Pipeline files (*.pipeline)|*.pipeline|All files (*.*)|*.*";
            openDialog.FilterIndex = 0;
            openDialog.RestoreDirectory = true;
            if (openDialog.ShowDialog() == DialogResult.OK)
            {
                FileStream stream = File.Open(openDialog.FileName, FileMode.Open);
                var formatter = new BinaryFormatter();
                while (stream.Position != stream.Length)
                {
                    var p = (Settings)formatter.Deserialize(stream);
                    Pipeline.Add(p);
                    PipNames.Add(p.name);
                    PipelineList.Items.Add(p.name);
                }
                stream.Close();
            }
        }

        private void FDRUpDown_ValueChanged(object sender, EventArgs e)
        {
            var str = FDRUpDown.Value.ToString();
            int digits = str.Substring(str.IndexOf(".")).Length;
            if (digits < 1) digits = 1;
            FDRUpDown.DecimalPlaces = digits;
        }

        private void PipUpdate_Click(object sender, EventArgs e)
        {
            if (PipelineList.SelectedIndices.Count == 1)
                if (!in_pipeline || PipelineList.SelectedIndices[0] > processed)
            {
                Settings S = new Settings();
                SaveSettings(ref S);
                Pipeline[PipelineList.SelectedIndices[0]] = S;
            }
        }

        private void PipelineList_MouseClick(object sender, MouseEventArgs e)
        {

        }

        private void PipelineList_SelectedIndexChanged(object sender, EventArgs e)
        {
            if (PipelineList.SelectedIndices.Count == 1)
            {
                LoadSettings(Pipeline[PipelineList.SelectedIndices[0]]);
            }
        }

        private void PipExec_Click(object sender, EventArgs e)
        {
            if (Pipeline.Count >= 1 && !running && !in_pipeline)
            {
                for (int i = 1; i < PipelineList.Items.Count; i++) PipelineList.Items[i].BackColor = System.Drawing.Color.White;
                processed = 0;
                in_pipeline = true;
                process.Dispose();
                RunProcess(false, Pipeline[0]);
            }
        }
    }
}
