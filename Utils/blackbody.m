
function cm = blackbody(n, opt_interp)
% Colormap: black body (perceptually uniform)
% Source:   https://www.kennethmoreland.com/color-advice/
% Converted for MatLab: Efe Ilicak based on code from Stéfan van der Walt,
% and Nathaniel Smith; table from Kenneth Moreland.

%-- Parse inputs ---------------------------------------------------------%
if ~exist('n', 'var'); n = []; end
if isempty(n)
   f = get(groot,'CurrentFigure');
   if isempty(f)
      n = size(get(groot,'DefaultFigureColormap'),1);
   else
      n = size(f.Colormap,1);
   end
end

% by default, interpolate in rgb space
if ~exist('opt_interp','var'); opt_interp = []; end
if isempty(opt_interp); opt_interp = 'rgb'; end
%-------------------------------------------------------------------------%


% data for colormap
cm = [
   0.013038855	0.003753703	0.002103028
0.02607771	0.007507407	0.004206056
0.039116565	0.01126111	0.006309084
0.051112818	0.015014814	0.008412112
0.061452019	0.018768517	0.01051514
0.070643263	0.02252222	0.012618168
0.078978066	0.026275924	0.014721196
0.086643619	0.030029627	0.016824224
0.093768357	0.03378333	0.018927251
0.100444796	0.037537034	0.021030279
0.106742123	0.041271634	0.023133307
0.112713635	0.044825979	0.025236335
0.118401394	0.048211428	0.027339363
0.123839262	0.051448139	0.029442391
0.129054968	0.054552615	0.031545419
0.134158812	0.05748457	0.033652416
0.139428264	0.060094809	0.035776007
0.144868037	0.062396588	0.037917415
0.150464634	0.064409436	0.040077014
0.156205677	0.066148416	0.042202264
0.162064287	0.067667256	0.044270886
0.167968117	0.06914956	0.046266761
0.173908795	0.070614116	0.048192562
0.179885077	0.072061075	0.050053053
0.185895824	0.07349058	0.051852438
0.191939989	0.07490276	0.053594455
0.198016608	0.076297738	0.055282434
0.204124791	0.077675625	0.056919365
0.210263715	0.079036524	0.058507936
0.216432613	0.080380529	0.060050578
0.222630776	0.081707728	0.061549494
0.228857539	0.0830182	0.06300669
0.235112281	0.084312017	0.064423993
0.241394423	0.085589245	0.065803079
0.247703418	0.086849943	0.067145483
0.254038752	0.088094165	0.068452618
0.260399943	0.089321957	0.069725783
0.266786533	0.090533361	0.070966182
0.27319809	0.091728413	0.072174925
0.279634203	0.092907144	0.073353045
0.286094483	0.094069579	0.074501496
0.292578561	0.095215739	0.075621169
0.299086082	0.09634564	0.076712892
0.305616711	0.097459293	0.077777437
0.312170126	0.098556705	0.078815524
0.318746019	0.099637877	0.079827827
0.325343194	0.100703071	0.080823047
0.331959604	0.101752811	0.081818152
0.338595148	0.102787073	0.08281315
0.345249722	0.103805825	0.08380805
0.351923219	0.104809036	0.084802862
0.358615535	0.10579667	0.085797595
0.36532656	0.106768688	0.086792257
0.372056188	0.10772505	0.087786857
0.378804309	0.10866571	0.088781404
0.385570814	0.109590621	0.089775905
0.392355595	0.110499732	0.090770369
0.399158541	0.111392991	0.091764803
0.405979545	0.112270341	0.092759215
0.412818496	0.113131721	0.093753612
0.419675287	0.113977069	0.094748002
0.426549809	0.114806319	0.095742392
0.433441956	0.115619402	0.096736788
0.440351618	0.116416245	0.097731198
0.447278691	0.117196773	0.098725628
0.454223067	0.117960906	0.099720084
0.461184643	0.118708562	0.100714573
0.468163313	0.119439656	0.101709101
0.475158975	0.120154097	0.102703674
0.482171524	0.120851792	0.103698298
0.489200859	0.121532646	0.104692979
0.49624688	0.122196557	0.105687723
0.503309485	0.122843421	0.106682534
0.510388575	0.123473131	0.107677418
0.517484053	0.124085573	0.108672381
0.52459582	0.124680633	0.109667428
0.531723779	0.125258189	0.110662564
0.538867836	0.125818116	0.111657794
0.546027896	0.126360287	0.112653122
0.553203864	0.126884567	0.113648554
0.560395648	0.127390817	0.114644094
0.567603156	0.127878896	0.115639747
0.574826298	0.128348655	0.116635517
0.582064982	0.12879994	0.117631408
0.589319121	0.129232594	0.118627426
0.596588625	0.129646454	0.119623573
0.603873408	0.130041349	0.120619854
0.611173384	0.130417104	0.121616274
0.618488467	0.13077354	0.122612836
0.625818572	0.131110468	0.123609544
0.633163617	0.131427695	0.124606401
0.640523518	0.131725022	0.125603412
0.647898194	0.13200224	0.126600581
0.655287563	0.132259136	0.12759791
0.662691547	0.132495488	0.128595403
0.670110065	0.132711067	0.129593064
0.677543039	0.132905637	0.130590896
0.684990393	0.13307895	0.131588902
0.692452049	0.133230754	0.132587086
0.699026944	0.135262857	0.133128782
0.702938503	0.142761643	0.132300872
0.706849161	0.150033679	0.131444721
0.710758934	0.157107264	0.130559597
0.71466784	0.164005743	0.129644729
0.718575896	0.170748624	0.12869931
0.722483119	0.177352385	0.127722488
0.726389524	0.183831082	0.126713367
0.730295128	0.190196808	0.125671003
0.734199946	0.196460049	0.124594398
0.738103995	0.202629958	0.123482497
0.742007289	0.208714575	0.122334186
0.745909844	0.214721004	0.121148281
0.749811675	0.220655552	0.119923525
0.753712796	0.226523847	0.118658584
0.757613224	0.232330931	0.117352036
0.76151297	0.238081341	0.116002363
0.765412051	0.243779172	0.114607944
0.769310481	0.249428135	0.113167043
0.773208272	0.255031605	0.111677797
0.777105438	0.260592655	0.110138203
0.781001994	0.266114098	0.108546102
0.784897953	0.271598509	0.106899164
0.788793327	0.277048256	0.105194862
0.79268813	0.28246552	0.103430456
0.796582375	0.287852311	0.101602958
0.800476073	0.29321049	0.099709107
0.804369239	0.298541781	0.097745327
0.808261883	0.303847783	0.095707686
0.812154019	0.309129985	0.09359184
0.816045659	0.314389772	0.091392974
0.819936813	0.319628435	0.089105724
0.823827495	0.324847182	0.086724086
0.827717715	0.330047143	0.084241303
0.831607485	0.335229374	0.081649725
0.835496816	0.340394868	0.078940635
0.83938572	0.345544556	0.076104027
0.843274208	0.350679315	0.073128323
0.84716229	0.355799968	0.070000002
0.851049977	0.360907292	0.066703114
0.854937281	0.366002019	0.06321862
0.858824211	0.371084842	0.059523482
0.862710777	0.376156412	0.055589391
0.866596991	0.381217349	0.051380902
0.870482862	0.386268237	0.046852647
0.8743684	0.39130963	0.041944975
0.878253616	0.396342055	0.036621943
0.882138518	0.401366009	0.031174102
0.886023118	0.406381966	0.025643215
0.889907423	0.411390377	0.020028782
0.891065382	0.418463568	0.021466471
0.891987123	0.425647055	0.02354603
0.892890899	0.432775622	0.025700955
0.893776626	0.439852242	0.027932126
0.894644216	0.446879675	0.030240426
0.895493582	0.45386049	0.032626735
0.896324635	0.46079708	0.035091935
0.897137287	0.467691677	0.037636908
0.897931446	0.474546369	0.040262536
0.898707022	0.481363111	0.042883775
0.899463921	0.488143734	0.045483033
0.90020205	0.49488996	0.048065003
0.900921314	0.501603405	0.050631303
0.901621618	0.508285591	0.053183372
0.902302864	0.514937951	0.055722496
0.902964954	0.52156184	0.058249831
0.903607789	0.528158534	0.060766417
0.904231269	0.534729241	0.063273193
0.904835291	0.541275104	0.065771011
0.905419752	0.547797206	0.068260644
0.905984549	0.554296573	0.070742798
0.906529575	0.560774178	0.073218117
0.907054723	0.567230948	0.075687191
0.907559886	0.57366776	0.078150561
0.908044952	0.580085452	0.080608723
0.908509812	0.586484819	0.083062138
0.908954352	0.59286662	0.085511229
0.909378458	0.599231579	0.087956386
0.909782015	0.605580384	0.090397972
0.910164905	0.611913696	0.092836324
0.91052701	0.618232142	0.095271754
0.910868209	0.624536325	0.097704553
0.91118838	0.63082682	0.100134991
0.911487399	0.637104177	0.102563322
0.911765142	0.643368923	0.104989782
0.91202148	0.649621563	0.10741459
0.912256285	0.655862581	0.109837955
0.912469425	0.662092442	0.112260069
0.91266077	0.668311592	0.114681114
0.912830183	0.674520457	0.11710126
0.912977528	0.680719449	0.119520669
0.913102667	0.686908962	0.12193949
0.913205459	0.693089376	0.124357866
0.913285762	0.699261055	0.126775931
0.91334343	0.705424352	0.12919381
0.913378318	0.711579603	0.131611623
0.913390276	0.717727134	0.134029482
0.913379152	0.72386726	0.136447493
0.913344793	0.730000281	0.138865757
0.913287044	0.736126489	0.141284369
0.913205744	0.742246166	0.143703418
0.913100735	0.748359581	0.14612299
0.912971851	0.754466997	0.148543166
0.912818928	0.760568665	0.150964022
0.912641796	0.766664831	0.153385631
0.912440284	0.772755728	0.155808062
0.912214218	0.778841584	0.15823138
0.911963421	0.78492262	0.160655647
0.911687714	0.790999049	0.163080923
0.911386913	0.797071074	0.165507263
0.911060832	0.803138896	0.167934722
0.910709283	0.809202707	0.170363349
0.910332073	0.815262693	0.172793194
0.909929007	0.821319035	0.175224302
0.909499887	0.827371906	0.177656716
0.909044509	0.833421477	0.18009048
0.908562669	0.839467911	0.182525632
0.908054157	0.845511367	0.18496221
0.907518761	0.851551999	0.187400251
0.906956262	0.857589956	0.189839789
0.906366442	0.863625383	0.192280856
0.905749074	0.86965842	0.194723484
0.905103931	0.875689203	0.197167703
0.90443078	0.881717865	0.199613542
0.903729384	0.887744534	0.202061026
0.9029995	0.893769334	0.204510183
0.902240883	0.899792386	0.206961036
0.905988894	0.904022129	0.234547643
0.912135367	0.907253839	0.272395932
0.918103004	0.910496798	0.307074851
0.923888987	0.913751437	0.339589276
0.929490474	0.917018177	0.370541693
0.934904589	0.92029743	0.400323761
0.94012842	0.923589601	0.429205126
0.945159007	0.926895085	0.45737956
0.949993335	0.930214267	0.484991026
0.954628329	0.933547527	0.512149353
0.959060844	0.936895233	0.538940162
0.963287659	0.940257747	0.565431413
0.967305469	0.943635421	0.591677866
0.971110875	0.947028602	0.617724221
0.974700378	0.950437625	0.643607371
0.978070368	0.953862819	0.669358057
0.981217117	0.957304506	0.695002116
0.98413677	0.960762999	0.720561419
0.986825332	0.964238604	0.746054603
0.989278661	0.967731617	0.771497643
0.991492452	0.971242331	0.7969043
0.993462231	0.974771028	0.822286483
0.995183339	0.978317983	0.847654544
0.996650916	0.981883466	0.873017509
0.997859893	0.985467736	0.898383277
0.998804968	0.989071049	0.923758779
0.999480596	0.992693652	0.949150115
0.999880966	0.996335784	0.974562663 
   ];


%-- Modify the colormap by interpolation ---------------------------------%
%   Note: Interpolation can be done in hsv or rgb space depending on opt_interp.
p = size(cm,1); % default size of colormap
if strcmp(opt_interp,'hsv')
    hsv = rgb2hsv(cm);
    hsv = interp1(1:p, hsv, linspace(1,p,n), 'linear');
    cm = hsv2rgb(hsv);
else
    cm = interp1(1:p, cm, linspace(1,p,n), 'linear');
end


end