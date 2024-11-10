use crate::quaternion::quadratic_order::{QuadraticOrder, QuadraticOrderEl};
use crate::util::Big;
use crate::{ec117, ec1757, ec214, ec509, util::bytes_from_str};

pub fn params32() {
    // p = 110985345351181940559839359037079551
    let A = ec117::Fq::new(
        &ec117::Fp::decode_reduce(&bytes_from_str("9827042539290988311061073437476016")),
        &ec117::Fp::decode_reduce(&bytes_from_str("74851534054032925543405982204421719")),
    );
    let EA = ec117::Curve::new(&A);

    let Px = ec117::Fq::new(
        &ec117::Fp::decode_reduce(&bytes_from_str("94208639506250591239460812981844796")),
        &ec117::Fp::decode_reduce(&bytes_from_str("2778322924810514545629313540375222")),
    );
    let Py = ec117::Fq::new(
        &ec117::Fp::decode_reduce(&bytes_from_str("52983347914014001368598403803362876")),
        &ec117::Fp::decode_reduce(&bytes_from_str("75608234133947709814825753280898001")),
    );

    let P = ec117::Point {
        X: Px,
        Y: Py,
        Z: ec117::Fq::ONE,
    };

    let Qx = ec117::Fq::new(
        &ec117::Fp::decode_reduce(&bytes_from_str("90055805557512327802376186357057367")),
        &ec117::Fp::decode_reduce(&bytes_from_str("109750291186252297193749404845454157")),
    );
    let Qy = ec117::Fq::new(
        &ec117::Fp::decode_reduce(&bytes_from_str("103764053627190847209017305966084403")),
        &ec117::Fp::decode_reduce(&bytes_from_str("93767318089512661131903081470994534")),
    );

    let Q = ec117::Point {
        X: Qx,
        Y: Qy,
        Z: ec117::Fq::ONE,
    };

    let omPx = ec117::Fq::new(
        &ec117::Fp::decode_reduce(&bytes_from_str("94656396045791334018841745121192250")),
        &ec117::Fp::decode_reduce(&bytes_from_str("38804522927790469747695470125130916")),
    );
    let omPy = ec117::Fq::new(
        &ec117::Fp::decode_reduce(&bytes_from_str("8616129072200153248015871509769894")),
        &ec117::Fp::decode_reduce(&bytes_from_str("40350091676899444441517261987606172")),
    );

    let omP = ec117::Point {
        X: omPx,
        Y: omPy,
        Z: ec117::Fq::ONE,
    };

    let omQx = ec117::Fq::new(
        &ec117::Fp::decode_reduce(&bytes_from_str("100870104592719437939350654020286413")),
        &ec117::Fp::decode_reduce(&bytes_from_str("26027362721692591519283159689746739")),
    );
    let omQy = ec117::Fq::new(
        &ec117::Fp::decode_reduce(&bytes_from_str("21326262706705111492955035179548229")),
        &ec117::Fp::decode_reduce(&bytes_from_str("53144258333902684409298433673474719")),
    );

    let omQ = ec117::Point {
        X: omQx,
        Y: omQy,
        Z: ec117::Fq::ONE,
    };

    let gen_denom = 2.big(); // quadratic order is generated by 1 and (1 + theta)/gen_denom.
    let field_disc = -11.big();
    let order_conductor = 19457.big();
    let quadratic_order =
        QuadraticOrder::new(field_disc.clone(), gen_denom, order_conductor.clone());

    // omega = 1/2*ϑ + 115517/2
    let omega = QuadraticOrderEl::new(115517.big(), 1.big(), 2.big(), quadratic_order.clone());
    // beta = 5/2*ϑ + 125787/2
    let beta = QuadraticOrderEl::new(125787.big(), 5.big(), 2.big(), quadratic_order.clone());

    let strategy: [usize; 103] = [
        48, 34, 21, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1,
        1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 14, 13, 8,
        5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 5, 3, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1,
    ];

    let two_dim = ec117::TwoDim::new(EA, omega, beta, P, Q, omP, omQ);
    let klapoti = ec117::Klapoti::new(quadratic_order, two_dim);

    let ideal = klapoti.secret();

    klapoti.act(ideal, 23, strategy.to_vec());
}

pub fn params64() {
    // p = 13575412597899949847778495756098141666107570891477033872626024447
    let A = ec214::Fq::new(
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "1292950946049871174637822016068695952313499707811732828713807401",
        )),
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "11638552349944361160217509745682746247479151631328891592249884628",
        )),
    );
    let EA = ec214::Curve::new(&A);

    let Px = ec214::Fq::new(
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "11723429863689620707436440690551082646129077314391124572415166271",
        )),
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "10456954105373228724166853278926964858170087793625218427799646343",
        )),
    );
    let Py = ec214::Fq::new(
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "6529457113899779665174483378227674168471749366522982109161668651",
        )),
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "2580221148911431123713889901557542302001805934065736789183452123",
        )),
    );

    let P = ec214::Point {
        X: Px,
        Y: Py,
        Z: ec214::Fq::ONE,
    };

    let Qx = ec214::Fq::new(
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "4846366763224230479763307074770335580637949267357089588784694136",
        )),
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "348055148369189936264026829801394744138471925808921710639879097",
        )),
    );
    let Qy = ec214::Fq::new(
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "9323569297698788639677461583991517933658585230294070071048851679",
        )),
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "2120762020776291485328773521530233377574257488787697652383584048",
        )),
    );

    let Q = ec214::Point {
        X: Qx,
        Y: Qy,
        Z: ec214::Fq::ONE,
    };

    let omPx = ec214::Fq::new(
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "1572040473296259226414846237835939480953848406791807972162041759",
        )),
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "6356558337637600514391792656026347990281988097418702304944537063",
        )),
    );
    let omPy = ec214::Fq::new(
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "3381572806759006732653322232692346340080480894439656638783798511",
        )),
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "5466589373957794996471797565399617192968835405114903341674905098",
        )),
    );

    let omP = ec214::Point {
        X: omPx,
        Y: omPy,
        Z: ec214::Fq::ONE,
    };

    let omQx = ec214::Fq::new(
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "4754960255176951360450680166150696950632914368595055674440383849",
        )),
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "11858848620189526315546198709763035166381181711595397756753332109",
        )),
    );
    let omQy = ec214::Fq::new(
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "5036852178299358060734894556022325626247115077702207703721822119",
        )),
        &ec214::Fp::decode_reduce(&bytes_from_str(
            "5132960215037716517189288106562229012491006181245242554710577441",
        )),
    );

    let omQ = ec214::Point {
        X: omQx,
        Y: omQy,
        Z: ec214::Fq::ONE,
    };

    let gen_denom = 2.big(); // quadratic order is generated by 1 and (1 + theta)/gen_denom.
    let field_disc = -11.big();
    let order_conductor = 469762049.big();
    let quadratic_order =
        QuadraticOrder::new(field_disc.clone(), gen_denom, order_conductor.clone());

    // omega = 1/2*ϑ + 1721613579/2
    let omega = QuadraticOrderEl::new(1721613579.big(), 1.big(), 2.big(), quadratic_order.clone());
    // beta = 7/2*ϑ + 21586103089/2
    let beta = QuadraticOrderEl::new(
        "21586103089".big(),
        7.big(),
        2.big(),
        quadratic_order.clone(),
    );

    let strategy: [usize; 202] = [
        79, 55, 34, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5,
        3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1,
        1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2,
        1, 1, 1, 1, 1, 24, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5,
        3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1,
        1, 8, 5, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1,
    ];

    let two_dim = ec214::TwoDim::new(EA, omega, beta, P, Q, omP, omQ);
    let klapoti = ec214::Klapoti::new(quadratic_order, two_dim);

    let ideal = klapoti.secret();

    klapoti.act(ideal, 39, strategy.to_vec());
}

pub fn params128() {
    // p = 864175120484581453683482079962486176185193500155369104423588921177379322250834082489183304374038697487834084609675858746433355728113743766078731283595263
    let A = ec509::Fq::new(
        &ec509::Fp::decode_reduce(&bytes_from_str("269704802029583887647805335061052403712442989875585143749461752225029821224697353454810243133823728633187546990793816614652026169443982629566764071566284")),
        &ec509::Fp::decode_reduce(&bytes_from_str("45075497449298242852915889845996934411298434793363088816673972940092720613863185721527461482928834564164321456798288718942619300938284615603658540248014")),
    );
    let EA = ec509::Curve::new(&A);

    let Px = ec509::Fq::new(
        &ec509::Fp::decode_reduce(&bytes_from_str("447892616027457465502233967364086561586562335469006266290465208783526124103064492970873055569941134579269272538807316962481005160382570106505154200728466")),
        &ec509::Fp::decode_reduce(&bytes_from_str("800290003970120280572593787002025355019721395023053350491644224698148995048995063007660746234392633217247174162428569404247459867752475229860556965198279")),
    );
    let Py = ec509::Fq::new(
        &ec509::Fp::decode_reduce(&bytes_from_str("519877388686727466920907141814943949020355617115239136230804736460919812194254261679112307050726660720836135386242885754976273452591631011423917422348771")),
        &ec509::Fp::decode_reduce(&bytes_from_str("605114052319867558932342087059807612022892998737365751462058895252259494372498773096167541799595985840585215351339563198397547806341378416496744576958128")),
    );

    let P = ec509::Point {
        X: Px,
        Y: Py,
        Z: ec509::Fq::ONE,
    };

    let Qx = ec509::Fq::new(
        &ec509::Fp::decode_reduce(&bytes_from_str("581507943031590973149523346726187302800330615020073945319870413336700595418472576727407722593565186779028791738389643633637621776018448806969800496171665")),
        &ec509::Fp::decode_reduce(&bytes_from_str("296992433797934585469880380198135514547268838524842360635380694206871872211981578239354201968494836926905151653846428203979911117308845627563297064024134")),
    );
    let Qy = ec509::Fq::new(
        &ec509::Fp::decode_reduce(&bytes_from_str("374104752400509605673980871207570386545470515033643307547656682897791939307273824348912788757046980260910734948072787715087017921252780682160868195437566")),
        &ec509::Fp::decode_reduce(&bytes_from_str("860886762997546709942891187910356315524716017603589228443181778475513536809307767636179123848468273049765258358389766247071465805595180262745557470648470")),
    );

    let Q = ec509::Point {
        X: Qx,
        Y: Qy,
        Z: ec509::Fq::ONE,
    };

    let omPx = ec509::Fq::new(
        &ec509::Fp::decode_reduce(&bytes_from_str("117364827899177059651277259150474935010104080982836495299275230154414640464707416846531424857872084495923618595737249168026438406721780781232254192922788")),
        &ec509::Fp::decode_reduce(&bytes_from_str("266058077597066160578639552635422053412206620907513149991022030262751902463737172585062035857693236118635861959941433228480542149175344816135702480159034")),
    );
    let omPy = ec509::Fq::new(
        &ec509::Fp::decode_reduce(&bytes_from_str("8022430588726828947198371836019131893534055152329420428850721991534462381340373664827471203284775511213176889211072082060440921574624796381053426098382")),
        &ec509::Fp::decode_reduce(&bytes_from_str("477238204785152495975270768173116194897686017354839920445991304698415818216147340967616942878422509035329189025930242592810529322598699646141543530614462")),
    );

    let omP = ec509::Point {
        X: omPx,
        Y: omPy,
        Z: ec509::Fq::ONE,
    };

    let omQx = ec509::Fq::new(
        &ec509::Fp::decode_reduce(&bytes_from_str("147775258981285866966116022022913564135666348961661433361136122574185027907525995802973409019591943944341051945210466497221185165048385369810117399652359")),
        &ec509::Fp::decode_reduce(&bytes_from_str("671827881927573715456912996097344332484172604565595244033691413943714532496526944794528845949791873076149479598009001017843295012167914989801350326034092")),
    );
    let omQy = ec509::Fq::new(
        &ec509::Fp::decode_reduce(&bytes_from_str("211785743333482542693275100259311165267101203331379805521488281141966213608336576673082845918500195026369271217051155063511774708462930773875347293706060")),
        &ec509::Fp::decode_reduce(&bytes_from_str("739118501506265494013422084991999169053986131949904071791478793677466824370959165027768700679545439036962047122800374196734886830129872014981877581011610")),
    );

    let omQ = ec509::Point {
        X: omQx,
        Y: omQy,
        Z: ec509::Fq::ONE,
    };

    let gen_denom = 2.big(); // quadratic order is generated by 1 and (1 + theta)/gen_denom.
    let field_disc = -11.big();
    let order_conductor = 469762049.big();
    let quadratic_order =
        QuadraticOrder::new(field_disc.clone(), gen_denom, order_conductor.clone());

    // omega = 1/2*ϑ + 3011552840461219095/2
    let omega = QuadraticOrderEl::new(
        "3011552840461219095".big(),
        1.big(),
        2.big(),
        quadratic_order.clone(),
    );
    // beta = 5/2*ϑ + 30751779181008197303/2
    let beta = QuadraticOrderEl::new(
        "30751779181008197303".big(),
        5.big(),
        2.big(),
        quadratic_order.clone(),
    );

    let strategy: [usize; 202] = [
        79, 55, 34, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5,
        3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1,
        1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2,
        1, 1, 1, 1, 1, 24, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5,
        3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1,
        1, 8, 5, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1,
    ];

    let two_dim = ec509::TwoDim::new(EA, omega, beta, P, Q, omP, omQ);
    let klapoti = ec509::Klapoti::new(quadratic_order, two_dim);

    let ideal = klapoti.secret();

    klapoti.act(ideal, 36, strategy.to_vec());
}

pub fn params512() {
    // p = 5076710506325114263887784071231455034756216907147249960253057351802083963386485171025944534774241392767171870098219516100819390082280320711738029767301300855033019711752666578955335646304863396566509330255028774784605158587986084272576435703062147836535217973939576405222935818910444049555313074695094377189706351825230979686667159010935188019400670487020900806597631804052117137991804103515765500501087078278636646850534420705487174308557136490766253617874076515911359437308290137101213848610281168195704356062374763729696849919
    let A = ec1757::Fq::new(
        &ec1757::Fp::decode_reduce(&bytes_from_str("1033497863051804585527435091442391243099210472035325657865292130796965301150496229942039090638314154092160983809839996833940571808214129411028962214661729077547573536119963581608635132947118932688378286113072662378221808696915431026273164775403679033037744761781351632380152860919397107675986945281917985945393038166260456910986734804444323729004203096789550173135834882120775946981779711437455654971503117578585095494044501886597156331547108178891795074234657755243861537588049782526471463970317563124421016348931558851564170003")),
        &ec1757::Fp::decode_reduce(&bytes_from_str("0")),
    );
    let EA = ec1757::Curve::new(&A);

    let Px = ec1757::Fq::new(
        &ec1757::Fp::decode_reduce(&bytes_from_str("5048910686393121072175511064161739598480795325508416182336961261480824537276601778369706715791149551876046150635228545567757343584949980128457097979274137675282069394224068034996091756267945689792845110496338052256320019713116359293729373068893066567556320583828301474535058539180536276063969004437635752827749304646515595260772759258525914170590213199633575048894206677433000754300469363222699602911832145887380787547136639868721389246926685938267279209733683181946165190702892774311392564619777849409241147789956967065930630855")),
        &ec1757::Fp::decode_reduce(&bytes_from_str("1848484917959088167713635648338125074974978488511506823357704532618801172411686040024569386768625133946910789364299776916675308309665171168731550637970362138621475471101150957966904717044395362566477458739214614295328306046704994204827651052966214998066576219878587757922991202416094486940375660753709171316079864405633649956211909242437412528826394171653894216156075864975204039744341712678512021542138583098949898528473644426237618196927481452405807087646048358381110914371375233117825256970953905774452722314061444263976138158")),
    );
    let Py = ec1757::Fq::new(
        &ec1757::Fp::decode_reduce(&bytes_from_str("1370618822146742373845778317881096270136754090099387370434661808560163313704265621424918429471559364333126498604896490297935793634965284877206995818915612124389342266963216337782780768212251316490516813059337452283455449208649349458882354542645429383116412334347877445180152978534954220558566966501610630412260191985214209853135035286204615366765705293934767314219433052211935629281471648289184199845202277251669170154286474034479632730721818425299461903284484730827536505704326767064054413282704219896800800039823078970406083681")),
        &ec1757::Fp::decode_reduce(&bytes_from_str("5021626645646545177459388742195734929500987090835331437981200130022835699587348557837066432024415992698431124728470850652731662510022368192422941638574789445959059551507810547689089418742674347835147035295051556124497545428213722291261073335642101405752027983719519805152008673468505834041292699385356234624162473691047708606239787373046235944618497814688293671515196982058625751177369413717201933585268592771035471087970418434444129058651637824810483043275701896269073070985577320205352425968624621321869556906178285575204614657")),
    );

    let P = ec1757::Point {
        X: Px,
        Y: Py,
        Z: ec1757::Fq::ONE,
    };

    let Qx = ec1757::Fq::new(
        &ec1757::Fp::decode_reduce(&bytes_from_str("779943298214243930133605267715764599036548565623401844149973989469087706305006822795138175697613669166419322414040396803328599265955906730707086145851912236773573333920282520239788064829358263962741752056085187833343776775772259808308158990411367049640303999594535219609813695778884644540196131979837739306841623791999062342007645846189568352818848672490094280750764441600740480689610400004692900619504969677132751673664559940436315913060307036401631682130994961845939857591219138536023973294466235281824419109094071504779841479")),
        &ec1757::Fp::decode_reduce(&bytes_from_str("2030385201556031938353848727773796221245835350773885731296974683677239293352689748638142958316973063602795176427636491504761683099884299236324018718432490788521466112732375356343715896840089786166769472392100009510717946728135850552465989697496656692174463579183495602269909974691656638332742579632074360531932705478316600424566730749876038027077190410663253993882215719757607024374206774766785863936171691450826832561759516491797854372926389571190267015820571103126883019918902038896462610949216866026274303781303512468035205989")),
    );
    let Qy = ec1757::Fq::new(
        &ec1757::Fp::decode_reduce(&bytes_from_str("3319410008078328993963229505907317976358418808487557417436970660838293498012211018768895993989884668698166269889094841842431649956209755518774150956118164758269783496265459195318291032444606440305769088081220143563694693417460431454597600988167613542789568792764904352080679718358607748339036798927914673481685045582093110560498276185366982818810139295876732258905417797679457809944522914477898526248920240100462800147682333925806708765980012631324692571742996115762204831579497438296681697873213263811018838757844519341876289129")),
        &ec1757::Fp::decode_reduce(&bytes_from_str("3600659067292362056902482530357307564484382189427643625642198634827624469348442900005933594138613625229571864229309921441909449853495364774707295329613061082061050527592840137506987734236779532670201530659376491076478963090223820189590902852126915727808691908676630452513224912502822623443118710829391820581969406361261921038472210958851674855312960113058900174755010975146736120432754521714767625528898671154330614896273076577615152727192996648579854335803358127441877432427169997493920467807366326624327190554861114805332413704")),
    );

    let Q = ec1757::Point {
        X: Qx,
        Y: Qy,
        Z: ec1757::Fq::ONE,
    };

    let omPx = ec1757::Fq::new(
        &ec1757::Fp::decode_reduce(&bytes_from_str("2529338630020879777912707502451458781199209414710588014331438074786516016581121018472944650290591691054408198532602897734239066800094362690333531601704545394794145786994747884837499208347539749195043864833545598775760570874920815872148804922531326608936120356021489328010662722740955041574718371954134474936793220033441736822674967436122682844116783828165402640222678234096184077374704289535067245233895039282669380407284242824959131725626988727698259121758031026683170882825042800744544888326510087786754494924461013478828735383")),
        &ec1757::Fp::decode_reduce(&bytes_from_str("2973284428789767884878777035721760522260893526395281734751721339051669226888714123868021925344704236387152634702650735805716311923601047415945284992337518732668510542061725725360486180734131552805462453490074039197350998406821755939516557841713304015805640468191212443144033639183412126590499204689647435853259951318935227616793978375530914659111328719226359810902251590438875109101158395005684293178819124620985374055225101366323085255905535162615206082583372637577283037111825748932879633995456069521196059974106362072839670295")),
    );
    let omPy = ec1757::Fq::new(
        &ec1757::Fp::decode_reduce(&bytes_from_str("3835261253847548761725798908755005498807823088085463056508345300366885607821803968839568548723465354138038602824234620537781652084144132864488344436157615297320088097147666832682099008475932038484033721870030272365298025459293893627711544604643748007705982943574361430922672310623769765086444729356212987412826783367302156522743831412758267384954504630709278960152303728977813260669631472380938295684025253936687946841881423274503606494441534686776210871386881473220465038322207986044743958560147927284670077151851056037492969873")),
        &ec1757::Fp::decode_reduce(&bytes_from_str("1417510832570480056088696074754544849467022966920840105779270617494414158671078130635062681845385240944653386659484028626131966167726200775251448507097636879975970302081871055142021969723284956480310739426241144174024163217220095741901344770858663857697355977248597256038145638393140218022573735537239540502566957072197804507976497420622569928094368472806368776461064997231720659143835675406988413561103932826846979535138934999449832987129043695925431817103441440901611631202067389977028565938384364059866776738876449477272303093")),
    );

    let omP = ec1757::Point {
        X: omPx,
        Y: omPy,
        Z: ec1757::Fq::ONE,
    };

    let omQx = ec1757::Fq::new(
        &ec1757::Fp::decode_reduce(&bytes_from_str("4692727323083168715187467510263579856800980163844460471870235745564104208841950029382106625947007772522950856020470440147064064191229179674182547967695764802764503419650272093955646602245380101101743155839880328465475785672069251388337716472161319054016303466449289876815405866921973117838245174891060289061220273801672622689176764586348664111945755096107996303450708413940386313739298021391995561263137155222924374528290190262560533097735564361523516699947370958507394886497929916440090420805055253474802658251620915984083606193")),
        &ec1757::Fp::decode_reduce(&bytes_from_str("3793949687840372072583152451655604075155490850752205943977747052868401379675720112965912188500563616687006933606899399272710199214874567919304412950771605331906198842183053403598558785258578333679691014339408666323700102297360548387979600205008844389176491993707521430633078477684367543073687657080121173956113385612778799434445520778818798287366954068952196382595667222917813696646223680486099202351549569663274627759600910724747319940179022141215748409193568601479189762364271524524127432895456688646668228488325874381173930003")),
    );
    let omQy = ec1757::Fq::new(
        &ec1757::Fp::decode_reduce(&bytes_from_str("2231972509438558837556988366057723035292451767991005729288654465395996359267013575431510092992657582601834258511746923600425138618858817103205054497403726373776369311411999224467605182831763364784728803895648826717471014106138895367493497612649889041850295646255616683098255855790421903385681365055168047694096503041781153649829956622042830679142799304634046179106282272681553372086820385434937222319231020737463459637135971284026626159061402895795648562533358383419314608065328306698076272188965961891633487740120169453744596383")),
        &ec1757::Fp::decode_reduce(&bytes_from_str("4781215979136101193408064689061727000548703137729904724150164933833931821519493341027055974238306097304894634804880366503712516355896196489384725655980550577530819549566554345138915738249311529217135435280112510993362781193223210894532020355220943175818643505765098397419263841734176656775384048577637941586971132875067499403293937031008389507638130266910823020614302132422789277468351144834828962684748656885216630397851056228859627202415116309556362248237798528786700656424155361933893306621717047577191047646113239931492735545")),
    );

    let omQ = ec1757::Point {
        X: omQx,
        Y: omQy,
        Z: ec1757::Fq::ONE,
    };

    let gen_denom = 2.big(); // quadratic order is generated by 1 and (1 + theta)/gen_denom.
    let field_disc = -11.big();
    let order_conductor =
        "34912628517411597986722280174925946878633333688546875202053332372871337315891".big();
    let quadratic_order =
        QuadraticOrder::new(field_disc.clone(), gen_denom, order_conductor.clone());

    // omega = 1/2*ϑ + 1329069277979427518548494869659323483948956584845223167970401259730478351270907/2
    let omega = QuadraticOrderEl::new(
        "1329069277979427518548494869659323483948956584845223167970401259730478351270907".big(),
        1.big(),
        2.big(),
        quadratic_order.clone(),
    );

    // beta = 27/2*ϑ + 1474947587453319227067450510036589234454539818780937894230669120639397875015835/2
    let beta = QuadraticOrderEl::new(
        "1474947587453319227067450510036589234454539818780937894230669120639397875015835".big(),
        27.big(),
        2.big(),
        quadratic_order.clone(),
    );

    let two_dim = ec1757::TwoDim::new(EA, omega, beta, P, Q, omP, omQ);
    let klapoti = ec1757::Klapoti::new(quadratic_order, two_dim);

    let ideal = klapoti.secret();

    let strategy: [usize; 1733] = [
        610, 377, 233, 157, 144, 89, 55, 34, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3,
        2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1,
        1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1,
        1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2,
        1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1,
        1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 55, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2,
        1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1,
        1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2,
        1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1,
        3, 2, 1, 1, 1, 1, 1, 55, 34, 21, 13, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1,
        1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3,
        2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1,
        1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1,
        1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 89, 55, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2,
        1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1,
        1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1,
        1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1,
        3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1,
        1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1,
        1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2,
        1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1,
        1, 1, 1, 1, 2, 1, 1, 1, 144, 89, 55, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3,
        2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5,
        3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1,
        1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2,
        1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5,
        3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1,
        1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1,
        2, 1, 1, 1, 55, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5,
        3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1,
        1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1,
        2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1,
        1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 233,
        144, 89, 55, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3,
        2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1,
        13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2,
        1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1,
        1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21,
        13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2,
        1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1,
        1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 55, 34, 21,
        13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2,
        1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1,
        1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5,
        3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1,
        8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 89, 55, 34, 21, 13, 8, 5, 3, 2,
        1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5,
        3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1,
        1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1,
        1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1,
        1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1,
        1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1,
        2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1,
        1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1,
    ];

    klapoti.act(ideal, 261, strategy.to_vec());
}

#[cfg(test)]
mod tests {
    use super::{params128, params32, params512, params64};
    use crate::util::core_cycles;

    #[test]
    fn params32_test() {
        params32();
    }

    #[test]
    fn params64_test() {
        params64();
    }

    #[test]
    fn params128_test() {
        params128();
    }

    #[test]
    fn params512_test() {
        params512();
    }

    #[test]
    fn cycles_test() {
        let mut tt = [0; 100];
        for i in 0..100 {
            let begin = core_cycles();
            params512();
            let end = core_cycles();
            tt[i] = end.wrapping_sub(begin);
        }
        tt.sort();
        println!(
            "Cycles:     {:?}",
            tt
        );
    }
}
