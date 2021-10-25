TCanvas* selection_plot(TTree* computed)
{
    TCanvas* c1 = new TCanvas();
    TH1F* all = new TH1F("all", "all", 113, 0.7, 12);
    TH1F* noflash = new TH1F("noflash", "noflash", 113, 0.7, 12);
    TH1F* nows = new TH1F("nows", "nows", 113, 0.7, 12);
    TH1F* prompt_withmuons = new TH1F("promptmu", "promptmu", 113, 0.7, 12);
    TH1F* delayed_withmuons = new TH1F("delayedmu", "delayedmu", 60, 6, 12);
    TH1F* prompts = new TH1F("prompts", "prompts", 113, 0.7, 12);
    TH1F* delayeds = new TH1F("delayeds", "delayeds", 60, 6, 12);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    all->SetLabelSize(0.05, "xyz");
    all->GetYaxis()->SetRangeUser(0.5, 2e7);
    all->GetYaxis()->SetTitle("Events");
    all->GetYaxis()->SetTitleSize(0.05);
    all->GetXaxis()->SetTitle("Energy (MeV)");
    all->GetXaxis()->SetTitleSize(0.05);
    int style = 20;
    int size = 1;

    all->SetMarkerStyle(style);
    all->SetMarkerColor(1);
    all->SetMarkerSize(size);

    noflash->SetMarkerStyle(style);
    noflash->SetMarkerColor(44);
    noflash->SetMarkerSize(size);

    nows->SetMarkerStyle(style);
    nows->SetMarkerColor(3);
    nows->SetMarkerSize(size);

    prompt_withmuons->SetMarkerStyle(style);
    prompt_withmuons->SetMarkerColor(4);
    prompt_withmuons->SetMarkerSize(size);

    delayed_withmuons->SetMarkerStyle(24);
    delayed_withmuons->SetMarkerColor(4);
    delayed_withmuons->SetMarkerSize(size);

    prompts->SetMarkerStyle(style);
    prompts->SetMarkerColor(2);
    prompts->SetMarkerSize(size);

    delayeds->SetMarkerStyle(24);
    delayeds->SetMarkerColor(2);
    delayeds->SetMarkerSize(size);

    computed->Draw("energy>>all", "detector == 1", "P");
    computed->Draw("energy >> noflash", "detector == 1 && (tag_flasher == 0 || tag_flasher == 2)", "same P");
    computed->Draw("energy >> nows", "detector == 1 && (tag_flasher == 0 || tag_flasher == 2) && !tag_WSMuonVeto", "same P");

    computed->Draw("energy_previous_PromptLike >> promptmu",
            "detector == 1 && !tag_flasher && !tag_WSMuonVeto && tag_DelayedLike && dt_previous_PromptLike < 200e3 && dt_previous_PromptLike > 1e3 && dt_next_DelayedLike > 200e3 && num_PromptLikes_400us == 1",
            "same P");

    computed->Draw("energy >> delayedmu", "detector == 1 && !tag_flasher && !tag_WSMuonVeto && tag_DelayedLike && dt_previous_PromptLike < 200e3 && dt_previous_PromptLike > 1e3 && dt_next_DelayedLike > 200e3 && num_PromptLikes_400us == 1", "same P");

    computed->Draw("energy_previous_PromptLike >> prompts", "detector == 1 && tag_IBDDelayed", "same P");
    computed->Draw("energy >> delayeds", "detector == 1 && tag_IBDDelayed", "same P");

    c1->SetLogy();

    c1->Modified();
    return c1;
}
