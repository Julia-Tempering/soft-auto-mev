#!/usr/bin/env -S julia --heap-size-hint=${task.memory.toGiga()}G
using Pkg
Pkg.activate(joinpath("$baseDir", "$julia_env")) 
include(joinpath("$baseDir", "$julia_env", "src", "utils.jl")) # loads dependencies too

@info """ Running experiments with the following combination of parameters:
	
	${arg}

"""

function main()
	# collect global vars 
	explorer_type = "${arg.sampler_type}"
	model = "${arg.model}"
    # collect pre tuned parameters
    if startswith(explorer_type, "RWMH")
        step_size = if model == "mRNA"
            0.0002012607744541562
        elseif model == "horseshoe_logit"
            0.016202325761796962
        elseif startswith(model, "funnel")
            0.4179156164274893
        elseif startswith(model, "banana")
            0.07696708866780329
        end
        estimated_target_std_deviations = if model == "mRNA"
            [0.0190223794994506, 0.026064173666313934, 0.03699569011340215, 0.017054173366615466, 0.02333280804373653]
        elseif model == "horseshoe_logit"
            [0.3976109369629807, 0.9984864862780988, 1.027494494666837, 0.830706966912766, 1.6652982655076982, 1.7008310729842921, 1.0196194879771403, 1.409644039311139, 1.8893578255615415, 1.0203390071458347, 1.4174477916907848, 1.2188698968889875, 0.9044355785762735, 1.1234600862598616, 1.650623112598901, 0.7241114503275085, 1.1489194093919959, 0.7995377704847795, 1.6954740775602086, 1.0687364428667059, 1.11701773303254, 1.1299705982285395, 0.7566719181072418, 1.055340558771287, 0.7308655277852215, 0.901036380388781, 0.9420660069406409, 0.907275903254109, 1.340487245223073, 0.946150226368304, 1.0634856282369578, 1.646652742775088, 0.5935021796009333, 0.8782853319180791, 0.8648404642387917, 0.8776177305015711, 1.0530501643616723, 0.7897111310037931, 0.9668070160372262, 1.0653676676364439, 1.5275319934481402, 1.6860589172862686, 1.3406045582194606, 1.971014090757069, 1.1943052557513028, 1.886860403167677, 1.489703709708733, 1.220952324276575, 1.013040914448645, 1.3856409385975985, 0.9643000483185463, 1.6742233671122155, 1.0533020313279793, 0.6854230501347263, 1.293683526951756, 0.711601122522722, 1.2975550946541832, 1.968665369446995, 1.329595405837782, 1.0896863693040093, 1.4575860394304563, 0.9352065644364521, 1.5126370763420127, 0.4227463597298082, 0.3556047974238951, 0.453095124822496, 1.8436195572170333, 0.9868895709880248, 1.0496084886870976, 0.9962514544222364, 1.363911152219638, 1.1318041445542197, 1.4374029040758929, 1.2606226966987195, 0.4373915370096597, 0.5009391115546739, 0.7568481741716551, 0.6352572164860639, 0.7361793741021262, 0.5723634647964634, 0.2198174995530747, 0.4458619831442608, 0.9695990649393453, 1.0717745088613246, 0.24222810089933106, 0.5301659169547552, 0.6989087086000637, 0.592710022221424, 0.39566525240265144, 0.45892903061645546, 0.5454124704914264, 0.9112005434381479, 0.5020637420540499, 0.8321316036018519, 0.576803120817221, 0.646028234777331, 0.48374834841750924, 0.6361658847872158, 0.7128342854664868, 0.5997977562546788, 0.6871939999695821, 1.3939655433439897, 0.3356997956321966, 0.8781848884774515, 1.8612390566382855, 1.5455702678615577, 1.8608625827609775, 3.1078963454184825, 0.47019894519168404, 3.5631355096491153, 0.3911254554006779, 0.3779042111328813, 0.4437305273828388, 1.0366604766993277, 1.1190760353905413, 0.6974043546489843, 0.4487929400214931, 0.439645509196368, 0.9970331439083179, 4.428272856463068, 0.27709609056623363, 0.7002408548126741]
        elseif startswith(model, "funnel")
            [2.7712018780501895, 116.8387973797307, 101.13224298644053, 89.71676443117684]
        elseif startswith(model, "banana")
            [1.8876761085451772, 6.166638145097786, 6.150015870527459, 6.145694093797878]
        end
    elseif startswith(explorer_type, "MALA")
        step_size = if model == "mRNA"
            0.030120243339653643
        elseif model == "horseshoe_logit"
            0.011675801015816842
        elseif startswith(model, "funnel")
            0.08331189981174276
        elseif startswith(model, "banana")
            0.10824917450735302
        end
        estimated_target_std_deviations = if model == "mRNA"
            [0.0580385819268161, 0.008516095252733325, 0.0181184233413692, 0.1681737852558746, 0.02535319686313943]
        elseif model == "horseshoe_logit"
            [0.2157391232048527, 0.7973676591161728, 1.0664901098106916, 0.573965461831839, 1.0769147772733083, 0.961776827870865, 0.831011053526207, 0.7605036574763175, 1.2696948044849392, 0.9328876199230531, 1.1168600838140852, 0.5413526422193677, 1.0213177227828099, 0.8349418567888971, 1.4703856347299549, 0.9702696917552582, 0.9120613423683323, 0.76396066089375, 0.9179464288128131, 0.9654936799690015, 0.7945328125514987, 0.6096173836484036, 0.6821891515354904, 0.9456950850352457, 1.169647376956703, 0.815932918004312, 0.8454362752452768, 0.9472477086675318, 0.9412778066613059, 0.9060538519344455, 0.6863323313533524, 0.8287770613297717, 0.6173802360858514, 0.8948902424904308, 1.0779805910866893, 0.6103003486272377, 0.7566494626896755, 1.058210125768933, 1.5020234144742308, 0.606189538397187, 0.9359324111294982, 1.1564680921562902, 0.68014885602684, 0.8063537593453531, 0.9968689923445706, 1.6263825079827332, 0.7196474483169015, 2.8134996006884694, 1.2758669351562173, 0.9700636167345679, 1.7837016486124704, 0.7917857968403158, 0.5895461225617256, 1.100849051018512, 0.6995105685548814, 0.4646902438253567, 0.7693798593594292, 0.7854009137703217, 0.8143427239996717, 0.8688347520190138, 0.7963123760884475, 0.7440334859864857, 1.3988394208487296, 0.25191698928695605, 0.5189670587474426, 6.284559298011592, 2.302477922259013, 1.1700550871586148, 1.0193834305426592, 0.9808805514954875, 0.6000972262607642, 1.134377474891778, 1.7252545876365324, 0.5178363981784497, 0.42428718751165, 0.5320126176942193, 0.9900098848866291, 0.6653721092689528, 1.1494831390843978, 0.9130629176121827, 0.8678882418168111, 0.7794698006925039, 0.7469526872314456, 0.6240898234595387, 1.2122585805424593, 1.089770543222063, 0.7206782158428372, 0.4347734358455377, 0.7501111481274326, 0.4859920670583844, 0.45973567534311144, 0.5842139589297597, 0.9384228838626689, 0.7314214439397767, 0.8098059476483971, 0.6816417736905354, 0.9211022754104848, 1.255124852875402, 1.0202623412600185, 0.6667417037139315, 0.7740139333812655, 1.380420040098556, 0.40133651761583294, 0.6078938206304726, 0.6014548870472816, 1.5801199176091343, 2.155798827454538, 2.749937699133312, 1.5045007927861074, 0.6161141810561043, 4.553447803311403, 1.1007859203410362, 1.4335894866871735, 0.9255180832548286, 1.5171722279721331, 0.6707561693912556, 0.6495792208903104, 1.5516670203052776, 1.482406955167004, 0.7391500611638822, 2.934904800072485, 0.8146032084575252]
        elseif startswith(model, "funnel")
            [2.059820729110829, 10.570466948716449, 15.191307221519663, 11.050446533982003]
        elseif startswith(model, "banana")
            [1.5937680905446259, 3.810635136051719, 3.815770626991041, 3.8138983231186367]
        end
    elseif startswith(explorer_type, "HMC")
        step_size = if model == "mRNA"
            0.0021292989352067557
        elseif model == "horseshoe_logit"
            0.019094114895321035
        elseif startswith(model, "funnel")
            0.03857081758647994
        elseif startswith(model, "banana")
            0.08738166108027126
        end
        estimated_target_std_deviations = if model == "mRNA"
            [0.03676740253376456, 0.004907361973147706, 0.006005114415121635, 0.24227547147234066, 0.024902748876005125]
        elseif model == "horseshoe_logit"
            [0.4262932832324372, 1.5546688925010068, 1.3781138384533493, 1.6166016053004397, 1.7831995220264172, 1.4736412201755755, 1.4258926015767504, 1.5322943204146884, 1.5744078388022644, 1.2203885177532463, 1.1839343025684235, 1.2844938066510652, 1.5651243540830468, 1.2554823790727658, 1.182345458445816, 1.4206697943553792, 1.2972908364216433, 1.2056969464763874, 1.2602661224185554, 1.075040378302276, 1.8310016341423938, 1.5059204515141527, 1.2694457472925198, 1.6292837916002145, 1.2041570503406767, 1.1311853784445902, 1.522401543688679, 1.2087433036714172, 1.1528082357086404, 1.109573826233193, 1.2799839848020933, 1.16636046952913, 1.3927576080394148, 1.2772038013620834, 1.0740178078023206, 1.4418169266678913, 1.2577943364564352, 1.4312175524376582, 1.3891856894937287, 1.32727175758735, 1.2348850816776644, 1.2196670836185854, 1.2246422833842963, 1.368117442595983, 1.3228716247749739, 1.2839660446267152, 1.5190462054850598, 1.510613706944478, 1.4611742464918969, 1.6622054323595103, 1.292198708552924, 1.3704750646924535, 1.8047178379860547, 1.1953399553755475, 1.3822823493447973, 1.4685620158992492, 1.6689779317441464, 1.9636888613410506, 1.4246929666842492, 1.2836899452510095, 1.698542222003043, 1.4143148136838763, 4.714950798624962, 1.8411613357608108, 2.1132660389841553, 7.294533087885655, 2.398136134738643, 2.021642427908645, 2.9358076435675304, 2.723521118067579, 1.47759457189973, 1.3362175885849152, 3.569803224288045, 1.956581577705165, 0.9775830531560613, 1.1401530848611205, 1.2094124422268437, 1.9985728327371506, 1.4001877907275897, 1.1746460175596292, 1.0212957805595702, 0.9233117700985243, 0.9535535911608999, 0.9162953539431316, 0.9915836317845155, 0.8329330781410791, 0.7405069927602039, 0.9048944644475905, 0.7406914820375099, 0.671217158503514, 0.714486491933021, 1.1567849713383809, 2.210601815323047, 1.2729785613770614, 1.0077407989486358, 1.0381049451520499, 1.0077328723872787, 1.9766333474108937, 1.8096624848032865, 1.2757193266709312, 1.218387007040236, 1.7804271659024962, 0.93762817930534, 0.9770286582784041, 1.5393872324113587, 2.546885971126121, 3.5223491355220142, 1.5025054045335613, 1.5181766718560712, 4.602808995144157, 9.485839348952412, 3.434402730866732, 2.531624704913411, 17.50988909646406, 1.33603354003735, 2.125908366386478, 2.8780662027055857, 5.330440743131294, 2.7930599992156546, 2.491830648209364, 1.9775933731894448, 1.7719926083030033]
        elseif startswith(model, "funnel")
            [2.3909040054192086, 9.082942470764872, 10.319316536542189, 10.349608331817024]
        elseif startswith(model, "banana")
            [1.9217285254281298, 5.030646686266274, 5.03782492724954, 5.0468280289327145]
        end
    end
    explorer = if explorer_type == "RWMH"
        SimpleRWMH(
            step_size = step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter())
        )
    elseif explorer_type == "MALA"
        SimpleAHMC(
            step_size = step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter()),
            int_time = AutoStep.FixedIntegrationTime()
        )
    elseif explorer_type == "HMC"
        SimpleAHMC(
            step_size = step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter()),
            int_time = AutoStep.AdaptiveFixedIntegrationTime()
        )
    elseif explorer_type == "RWMH0.1"
        SimpleRWMH(
            step_size = 0.1 * step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter())
        )
    elseif explorer_type == "MALA0.1"
        SimpleAHMC(
            step_size = 0.1 * step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter()),
            int_time = AutoStep.FixedIntegrationTime()
        )
    elseif explorer_type == "HMC0.1"
        SimpleAHMC(
            step_size = 0.1 * step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter()),
            int_time = AutoStep.AdaptiveFixedIntegrationTime()
        )
    elseif explorer_type == "RWMH10.0"
        SimpleRWMH(
            step_size = 10.0 * step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter())
        )
    elseif explorer_type == "MALA10.0"
        SimpleAHMC(
            step_size = 10.0 * step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter()),
            int_time = AutoStep.FixedIntegrationTime()
        )
    elseif explorer_type == "HMC10.0"
        SimpleAHMC(
            step_size = 10.0 * step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter()),
            int_time = AutoStep.AdaptiveFixedIntegrationTime()
        )
    elseif explorer_type == "RWMH4.0"
        SimpleRWMH(
            step_size = 4.0 * step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter())
        )
    elseif explorer_type == "MALA4.0"
        SimpleAHMC(
            step_size = 4.0 * step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter()),
            int_time = AutoStep.FixedIntegrationTime()
        )
    elseif explorer_type == "HMC4.0"
        SimpleAHMC(
            step_size = 4.0 * step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter()),
            int_time = AutoStep.AdaptiveFixedIntegrationTime()
        )
    elseif explorer_type == "RWMH0.25"
        SimpleRWMH(
            step_size = 0.25 * step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter())
        )
    elseif explorer_type == "MALA0.25"
        SimpleAHMC(
            step_size = 0.25 * step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter()),
            int_time = AutoStep.FixedIntegrationTime()
        )
    elseif explorer_type == "HMC0.25"
        SimpleAHMC(
            step_size = 0.25 * step_size, 
            estimated_target_std_deviations = estimated_target_std_deviations,
            step_size_selector = AutoStep.ASNonAdaptiveSelector(),
            step_jitter = AutoStep.StepJitter(dist = Dirac(0.0), adapt_strategy = AutoStep.FixedStepJitter()),
            int_time = AutoStep.AdaptiveFixedIntegrationTime()
        )
    else
        make_explorer(explorer_type, "${arg.selector}", "${arg.int_time}", "${arg.logstep_jitter}")
    end
    target = if model == "funnel(4,1)"
            Pigeons.stan_funnel(3, 1.0)
        elseif model == "banana(4,1)"
            Pigeons.stan_banana(3, 1.0)
        elseif model == "funnel(4,0.3)"
            Pigeons.stan_funnel(3, 0.3)
        elseif model == "banana(4,0.3)"
            Pigeons.stan_banana(3, 0.3)
        else
            startswith(model, "horseshoe") ? stan_logpotential(model; dataset = "sonar") : stan_logpotential(model)
        end
	seed = ${arg.seed}
	miness_threshold = ${params.dryRun ? 1 : 100}

	samples, stats_df = if miness_threshold == 1 
        pt_sample_from_model(model, target, seed, explorer, miness_threshold)
    else
        pt_sample_from_model_fixed(model, target, seed, explorer, 19)
    end

	isdir("csvs") || mkdir("csvs")
	CSV.write("csvs/summary.csv", stats_df)
end

function Pigeons.initialization(target::StanLogPotential{StanModel,String,Pigeons.Immutable{String},Nothing}, 
    rng::AbstractRNG, ::Int64)
    d_unc = BridgeStan.param_unc_num(target.model) # number of unconstrained parameters
    init = randn(rng, d_unc)
    return Pigeons.StanState(init, StanRNG(target.model, rand(rng, UInt32)))
end

main()
