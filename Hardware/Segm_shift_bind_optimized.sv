module skeleton_final_optim1 # (
    parameter integer D = 1024,
    //parameter real p_sparse = 0.0078125, //0.0078125,
    parameter integer E = 64,
    parameter integer NB_SEGMENTS = 8,
    parameter integer LENGTH_SEGMENT = 128,
    parameter integer NB_TO_BUNDLE_IN_TIME = 256,
    parameter integer LBP_LENGTH = 6,
    parameter integer NB_CLASSES = 2,
    parameter integer vector_fold_factor = 1, //either 1,2,4...NB_SEGMENTS
    parameter integer channel_fold_factor = {CFF} //either 1,2,4...E (if E is a power of 2)
)
(
    input  clk_i,
    input  arst_n_in,
    input  [LBP_LENGTH-1:0] LBP_codes63,
    input  [LBP_LENGTH-1:0] LBP_codes62,
    input  [LBP_LENGTH-1:0] LBP_codes61,
    input  [LBP_LENGTH-1:0] LBP_codes60,
    input  [LBP_LENGTH-1:0] LBP_codes59,
    input  [LBP_LENGTH-1:0] LBP_codes58,
    input  [LBP_LENGTH-1:0] LBP_codes57,
    input  [LBP_LENGTH-1:0] LBP_codes56,
    input  [LBP_LENGTH-1:0] LBP_codes55,
    input  [LBP_LENGTH-1:0] LBP_codes54,
    input  [LBP_LENGTH-1:0] LBP_codes53,
    input  [LBP_LENGTH-1:0] LBP_codes52,
    input  [LBP_LENGTH-1:0] LBP_codes51,
    input  [LBP_LENGTH-1:0] LBP_codes50,
    input  [LBP_LENGTH-1:0] LBP_codes49,
    input  [LBP_LENGTH-1:0] LBP_codes48,
    input  [LBP_LENGTH-1:0] LBP_codes47,
    input  [LBP_LENGTH-1:0] LBP_codes46,
    input  [LBP_LENGTH-1:0] LBP_codes45,
    input  [LBP_LENGTH-1:0] LBP_codes44,
    input  [LBP_LENGTH-1:0] LBP_codes43,
    input  [LBP_LENGTH-1:0] LBP_codes42,
    input  [LBP_LENGTH-1:0] LBP_codes41,
    input  [LBP_LENGTH-1:0] LBP_codes40,
    input  [LBP_LENGTH-1:0] LBP_codes39,
    input  [LBP_LENGTH-1:0] LBP_codes38,
    input  [LBP_LENGTH-1:0] LBP_codes37,
    input  [LBP_LENGTH-1:0] LBP_codes36,
    input  [LBP_LENGTH-1:0] LBP_codes35,
    input  [LBP_LENGTH-1:0] LBP_codes34,
    input  [LBP_LENGTH-1:0] LBP_codes33,
    input  [LBP_LENGTH-1:0] LBP_codes32,
    input  [LBP_LENGTH-1:0] LBP_codes31,
    input  [LBP_LENGTH-1:0] LBP_codes30,
    input  [LBP_LENGTH-1:0] LBP_codes29,
    input  [LBP_LENGTH-1:0] LBP_codes28,
    input  [LBP_LENGTH-1:0] LBP_codes27,
    input  [LBP_LENGTH-1:0] LBP_codes26,
    input  [LBP_LENGTH-1:0] LBP_codes25,
    input  [LBP_LENGTH-1:0] LBP_codes24,
    input  [LBP_LENGTH-1:0] LBP_codes23,
    input  [LBP_LENGTH-1:0] LBP_codes22,
    input  [LBP_LENGTH-1:0] LBP_codes21,
    input  [LBP_LENGTH-1:0] LBP_codes20,
    input  [LBP_LENGTH-1:0] LBP_codes19,
    input  [LBP_LENGTH-1:0] LBP_codes18,
    input  [LBP_LENGTH-1:0] LBP_codes17,
    input  [LBP_LENGTH-1:0] LBP_codes16,
    input  [LBP_LENGTH-1:0] LBP_codes15,
    input  [LBP_LENGTH-1:0] LBP_codes14,
    input  [LBP_LENGTH-1:0] LBP_codes13,
    input  [LBP_LENGTH-1:0] LBP_codes12,
    input  [LBP_LENGTH-1:0] LBP_codes11,
    input  [LBP_LENGTH-1:0] LBP_codes10,
    input  [LBP_LENGTH-1:0] LBP_codes9,
    input  [LBP_LENGTH-1:0] LBP_codes8,
    input  [LBP_LENGTH-1:0] LBP_codes7,
    input  [LBP_LENGTH-1:0] LBP_codes6,
    input  [LBP_LENGTH-1:0] LBP_codes5,
    input  [LBP_LENGTH-1:0] LBP_codes4,
    input  [LBP_LENGTH-1:0] LBP_codes3,
    input  [LBP_LENGTH-1:0] LBP_codes2,
    input  [LBP_LENGTH-1:0] LBP_codes1,
    input  [LBP_LENGTH-1:0] LBP_codes0,
    input  [0:D-1] ictal_hv_in,
    input  [0:D-1] interictal_hv_in,
    output classification_ready,
    output [$clog2(NB_CLASSES)-1:0] classification,
    output send_next_LBP,
    output [$clog2(D):0] ictal_sim,
    output [$clog2(D):0] interictal_sim
);

//LBP_codes
logic [LBP_LENGTH-1:0] LBP_codes [E-1:0];
assign LBP_codes[63] = LBP_codes63;
assign LBP_codes[62] = LBP_codes62;
assign LBP_codes[61] = LBP_codes61;
assign LBP_codes[60] = LBP_codes60;
assign LBP_codes[59] = LBP_codes59;
assign LBP_codes[58] = LBP_codes58;
assign LBP_codes[57] = LBP_codes57;
assign LBP_codes[56] = LBP_codes56;
assign LBP_codes[55] = LBP_codes55;
assign LBP_codes[54] = LBP_codes54;
assign LBP_codes[53] = LBP_codes53;
assign LBP_codes[52] = LBP_codes52;
assign LBP_codes[51] = LBP_codes51;
assign LBP_codes[50] = LBP_codes50;
assign LBP_codes[49] = LBP_codes49;
assign LBP_codes[48] = LBP_codes48;
assign LBP_codes[47] = LBP_codes47;
assign LBP_codes[46] = LBP_codes46;
assign LBP_codes[45] = LBP_codes45;
assign LBP_codes[44] = LBP_codes44;
assign LBP_codes[43] = LBP_codes43;
assign LBP_codes[42] = LBP_codes42;
assign LBP_codes[41] = LBP_codes41;
assign LBP_codes[40] = LBP_codes40;
assign LBP_codes[39] = LBP_codes39;
assign LBP_codes[38] = LBP_codes38;
assign LBP_codes[37] = LBP_codes37;
assign LBP_codes[36] = LBP_codes36;
assign LBP_codes[35] = LBP_codes35;
assign LBP_codes[34] = LBP_codes34;
assign LBP_codes[33] = LBP_codes33;
assign LBP_codes[32] = LBP_codes32;
assign LBP_codes[31] = LBP_codes31;
assign LBP_codes[30] = LBP_codes30;
assign LBP_codes[29] = LBP_codes29;
assign LBP_codes[28] = LBP_codes28;
assign LBP_codes[27] = LBP_codes27;
assign LBP_codes[26] = LBP_codes26;
assign LBP_codes[25] = LBP_codes25;
assign LBP_codes[24] = LBP_codes24;
assign LBP_codes[23] = LBP_codes23;
assign LBP_codes[22] = LBP_codes22;
assign LBP_codes[21] = LBP_codes21;
assign LBP_codes[20] = LBP_codes20;
assign LBP_codes[19] = LBP_codes19;
assign LBP_codes[18] = LBP_codes18;
assign LBP_codes[17] = LBP_codes17;
assign LBP_codes[16] = LBP_codes16;
assign LBP_codes[15] = LBP_codes15;
assign LBP_codes[14] = LBP_codes14;
assign LBP_codes[13] = LBP_codes13;
assign LBP_codes[12] = LBP_codes12;
assign LBP_codes[11] = LBP_codes11;
assign LBP_codes[10] = LBP_codes10;
assign LBP_codes[9] = LBP_codes9;
assign LBP_codes[8] = LBP_codes8;
assign LBP_codes[7] = LBP_codes7;
assign LBP_codes[6] = LBP_codes6;
assign LBP_codes[5] = LBP_codes5;
assign LBP_codes[4] = LBP_codes4;
assign LBP_codes[3] = LBP_codes3;
assign LBP_codes[2] = LBP_codes2;
assign LBP_codes[1] = LBP_codes1;
assign LBP_codes[0] = LBP_codes0;



logic [$clog2(vector_fold_factor)-1:0] vf_counter;
logic [$clog2(channel_fold_factor)-1:0] cf_counter;

//counters for vector folding and channel folding
always @(posedge clk_i) begin
    if ((arst_n_in == 0) | ((vector_fold_factor == 1) & (channel_fold_factor == 1))) begin
        vf_counter <= 'd0;
        cf_counter <= 'd0;
    end
    else if (vector_fold_factor == 1) begin
        vf_counter <= 'd0;
        cf_counter <= cf_counter+1;
    end
    else if (channel_fold_factor == 1) begin
        vf_counter <= (cf_counter == (channel_fold_factor-1)) ? vf_counter+1:vf_counter;
        cf_counter <= 'd0;
    end
    else begin
        vf_counter <= (cf_counter == (channel_fold_factor-1)) ? vf_counter+1:vf_counter;
        cf_counter <= cf_counter+1;
    end
end

//signal that need next LBP_codes
assign send_next_LBP = ((vf_counter == (vector_fold_factor-1)) & (cf_counter == (channel_fold_factor-1))) ? 1:0;

logic [LBP_LENGTH-1:0] LBP_codes_saved [E-1:0];
logic [LBP_LENGTH-1:0] LBP_codes_current [E-1:0];
//regs for LBP_values, needed if channel or vector folding is used
int i;
always_comb begin
    if ((vf_counter > 0) | (cf_counter > 0)) begin
        for (i=0;i<E;i=i+1) begin
            LBP_codes_current[i] = LBP_codes_saved[i];
        end
    end
    else begin
        for (i=0;i<E;i=i+1) begin
            LBP_codes_current[i] = LBP_codes[i];
        end
    end
end

int i2;
always @(posedge clk_i) begin
for (i2=0;i2<E;i2=i2+1) begin
    if (arst_n_in == 0) begin
        LBP_codes_saved[i2] <= 'd0;
    end
    else begin
        LBP_codes_saved[i2] <= ((vf_counter == 0) & (cf_counter == 0)) ? LBP_codes[i2]:LBP_codes_saved[i2];
    end
end
end


logic [E/channel_fold_factor-1:0][(NB_SEGMENTS/vector_fold_factor)-1:0][$clog2(LENGTH_SEGMENT)-1:0] LBP_hvs;

logic [E/channel_fold_factor-1:0][0:LENGTH_SEGMENT*NB_SEGMENTS/vector_fold_factor-1] bind_outputs;

logic bundled_hv [0:LENGTH_SEGMENT*NB_SEGMENTS/vector_fold_factor-1];

logic [0:D-1] result_hv;
logic [$clog2(NB_TO_BUNDLE_IN_TIME):0] counter_bund_time;
logic [$clog2(NB_TO_BUNDLE_IN_TIME)-1:0] counter_array [0:D-1];

logic [0:D-1] ictal_hv;
logic [0:D-1] interictal_hv;


assign ictal_hv = (arst_n_in) ? ictal_hv : ictal_hv_in;
assign interictal_hv = (arst_n_in) ? interictal_hv : interictal_hv_in;

//The Big Five
IM_o1 #(.D(D), /*.p_sparse(p_sparse),*/ .E(E), .NB_SEGMENTS(NB_SEGMENTS), .LENGTH_SEGMENT(LENGTH_SEGMENT), .LBP_LENGTH(LBP_LENGTH), .vector_fold_factor(vector_fold_factor), .channel_fold_factor(channel_fold_factor)) 
IM_module (.clk_i(clk_i), .arst_n_in(arst_n_in), .LBP_codes(LBP_codes_current), .vf_counter(vf_counter), .cf_counter(cf_counter), .LBP_hvs(LBP_hvs));

BINDING_o1 #(.D(D), /*.p_sparse(p_sparse),*/ .E(E), .NB_SEGMENTS(NB_SEGMENTS), .LENGTH_SEGMENT(LENGTH_SEGMENT), .LBP_LENGTH(LBP_LENGTH), .vector_fold_factor(vector_fold_factor), .channel_fold_factor(channel_fold_factor)) 
BINDING_module (.LBP_hvs(LBP_hvs), .vf_counter(vf_counter), .cf_counter(cf_counter), .bind_outputs(bind_outputs));

BUNDLING_SPACE_o1 #(.D(D), /*.p_sparse(p_sparse),*/ .E(E), .NB_SEGMENTS(NB_SEGMENTS), .LENGTH_SEGMENT(LENGTH_SEGMENT), .LBP_LENGTH(LBP_LENGTH), .vector_fold_factor(vector_fold_factor), .channel_fold_factor(channel_fold_factor)) 
BUNDLING_SPACE_module (.clk_i(clk_i), .arst_n_in(arst_n_in), .bind_outputs(bind_outputs), .cf_counter(cf_counter), .bundled_hv(bundled_hv));

BUNDLING_TIME_o1 #(.D(D), /*.p_sparse(p_sparse),*/ .E(E), .NB_SEGMENTS(NB_SEGMENTS), .LENGTH_SEGMENT(LENGTH_SEGMENT), .NB_TO_BUNDLE_IN_TIME(NB_TO_BUNDLE_IN_TIME), .LBP_LENGTH(LBP_LENGTH), .vector_fold_factor(vector_fold_factor), .channel_fold_factor(channel_fold_factor)) 
BUNDLING_TIME_module (.clk_i(clk_i), .arst_n_in(arst_n_in), .bundled_hv(bundled_hv), .vf_counter(vf_counter), .cf_counter(cf_counter), .result_hv(result_hv), .counter_bund_time(counter_bund_time), .counter_array(counter_array));

SIMILARITY_SEARCH_o1 #(.D(D), /*.p_sparse(p_sparse),*/ .E(E), .NB_SEGMENTS(NB_SEGMENTS), .LENGTH_SEGMENT(LENGTH_SEGMENT), .NB_TO_BUNDLE_IN_TIME(NB_TO_BUNDLE_IN_TIME), .LBP_LENGTH(LBP_LENGTH), .vector_fold_factor(vector_fold_factor), .channel_fold_factor(channel_fold_factor), .NB_CLASSES(NB_CLASSES)) 
SIMILARITY_SEARCH_module (.clk_i(clk_i), .arst_n_in(arst_n_in), .result_hv(result_hv), .ictal_hv(ictal_hv), .interictal_hv(interictal_hv), .vf_counter(vf_counter), .cf_counter(cf_counter), .counter_bund_time(counter_bund_time), .classification(classification), .classification_ready(classification_ready), .ictal_sim(ictal_sim), .interictal_sim(interictal_sim));

endmodule


module IM_o1 # (
    parameter integer D = 1024,
    //parameter real p_sparse = 0.0078125,
    parameter integer E = 64,
    parameter integer NB_SEGMENTS = 8,
    parameter integer LENGTH_SEGMENT = 128,
    parameter integer LBP_LENGTH = 6,
    parameter integer vector_fold_factor = 1, //either 1,2,4...NB_SEGMENTS
    parameter integer channel_fold_factor = 1 //either 1,2,4...E (if E is a power of 2)
)
(
    input  clk_i,
    input  arst_n_in,
    input  [LBP_LENGTH-1:0] LBP_codes [E-1:0],
    input  [$clog2(vector_fold_factor)-1:0] vf_counter,
    input  [$clog2(channel_fold_factor)-1:0] cf_counter,
    output logic [E/channel_fold_factor-1:0][(NB_SEGMENTS/vector_fold_factor)-1:0][$clog2(LENGTH_SEGMENT)-1:0] LBP_hvs
);

//IM
logic [2**(LBP_LENGTH)-1:0][NB_SEGMENTS-1:0][$clog2(LENGTH_SEGMENT)-1:0] IM;

assign IM[0][0] = 'd1; assign IM[0][1] = 'd75; assign IM[0][2] = 'd60; assign IM[0][3] = 'd76; assign IM[0][4] = 'd85; assign IM[0][5] = 'd63; assign IM[0][6] = 'd104; assign IM[0][7] = 'd73; 
assign IM[1][0] = 'd56; assign IM[1][1] = 'd103; assign IM[1][2] = 'd117; assign IM[1][3] = 'd12; assign IM[1][4] = 'd25; assign IM[1][5] = 'd74; assign IM[1][6] = 'd69; assign IM[1][7] = 'd96; 
assign IM[2][0] = 'd90; assign IM[2][1] = 'd116; assign IM[2][2] = 'd50; assign IM[2][3] = 'd1; assign IM[2][4] = 'd27; assign IM[2][5] = 'd30; assign IM[2][6] = 'd54; assign IM[2][7] = 'd42; 
assign IM[3][0] = 'd44; assign IM[3][1] = 'd1; assign IM[3][2] = 'd23; assign IM[3][3] = 'd99; assign IM[3][4] = 'd117; assign IM[3][5] = 'd109; assign IM[3][6] = 'd120; assign IM[3][7] = 'd29; 
assign IM[4][0] = 'd40; assign IM[4][1] = 'd18; assign IM[4][2] = 'd67; assign IM[4][3] = 'd73; assign IM[4][4] = 'd36; assign IM[4][5] = 'd75; assign IM[4][6] = 'd1; assign IM[4][7] = 'd122; 
assign IM[5][0] = 'd127; assign IM[5][1] = 'd85; assign IM[5][2] = 'd5; assign IM[5][3] = 'd123; assign IM[5][4] = 'd89; assign IM[5][5] = 'd64; assign IM[5][6] = 'd59; assign IM[5][7] = 'd53; 
assign IM[6][0] = 'd64; assign IM[6][1] = 'd24; assign IM[6][2] = 'd117; assign IM[6][3] = 'd117; assign IM[6][4] = 'd30; assign IM[6][5] = 'd100; assign IM[6][6] = 'd53; assign IM[6][7] = 'd97; 
assign IM[7][0] = 'd3; assign IM[7][1] = 'd6; assign IM[7][2] = 'd47; assign IM[7][3] = 'd61; assign IM[7][4] = 'd48; assign IM[7][5] = 'd119; assign IM[7][6] = 'd80; assign IM[7][7] = 'd97; 
assign IM[8][0] = 'd44; assign IM[8][1] = 'd0; assign IM[8][2] = 'd72; assign IM[8][3] = 'd124; assign IM[8][4] = 'd39; assign IM[8][5] = 'd82; assign IM[8][6] = 'd41; assign IM[8][7] = 'd62; 
assign IM[9][0] = 'd70; assign IM[9][1] = 'd29; assign IM[9][2] = 'd116; assign IM[9][3] = 'd12; assign IM[9][4] = 'd37; assign IM[9][5] = 'd85; assign IM[9][6] = 'd61; assign IM[9][7] = 'd60; 
assign IM[10][0] = 'd109; assign IM[10][1] = 'd56; assign IM[10][2] = 'd43; assign IM[10][3] = 'd113; assign IM[10][4] = 'd105; assign IM[10][5] = 'd126; assign IM[10][6] = 'd3; assign IM[10][7] = 'd112; 
assign IM[11][0] = 'd91; assign IM[11][1] = 'd5; assign IM[11][2] = 'd124; assign IM[11][3] = 'd34; assign IM[11][4] = 'd102; assign IM[11][5] = 'd81; assign IM[11][6] = 'd123; assign IM[11][7] = 'd48; 
assign IM[12][0] = 'd115; assign IM[12][1] = 'd67; assign IM[12][2] = 'd49; assign IM[12][3] = 'd58; assign IM[12][4] = 'd115; assign IM[12][5] = 'd85; assign IM[12][6] = 'd70; assign IM[12][7] = 'd97; 
assign IM[13][0] = 'd94; assign IM[13][1] = 'd16; assign IM[13][2] = 'd111; assign IM[13][3] = 'd119; assign IM[13][4] = 'd48; assign IM[13][5] = 'd108; assign IM[13][6] = 'd30; assign IM[13][7] = 'd47; 
assign IM[14][0] = 'd120; assign IM[14][1] = 'd88; assign IM[14][2] = 'd116; assign IM[14][3] = 'd31; assign IM[14][4] = 'd2; assign IM[14][5] = 'd10; assign IM[14][6] = 'd44; assign IM[14][7] = 'd11; 
assign IM[15][0] = 'd107; assign IM[15][1] = 'd109; assign IM[15][2] = 'd80; assign IM[15][3] = 'd68; assign IM[15][4] = 'd73; assign IM[15][5] = 'd24; assign IM[15][6] = 'd111; assign IM[15][7] = 'd124; 
assign IM[16][0] = 'd12; assign IM[16][1] = 'd50; assign IM[16][2] = 'd118; assign IM[16][3] = 'd39; assign IM[16][4] = 'd117; assign IM[16][5] = 'd36; assign IM[16][6] = 'd2; assign IM[16][7] = 'd98; 
assign IM[17][0] = 'd120; assign IM[17][1] = 'd83; assign IM[17][2] = 'd51; assign IM[17][3] = 'd88; assign IM[17][4] = 'd106; assign IM[17][5] = 'd113; assign IM[17][6] = 'd57; assign IM[17][7] = 'd40; 
assign IM[18][0] = 'd7; assign IM[18][1] = 'd1; assign IM[18][2] = 'd123; assign IM[18][3] = 'd6; assign IM[18][4] = 'd23; assign IM[18][5] = 'd21; assign IM[18][6] = 'd19; assign IM[18][7] = 'd88; 
assign IM[19][0] = 'd40; assign IM[19][1] = 'd33; assign IM[19][2] = 'd98; assign IM[19][3] = 'd36; assign IM[19][4] = 'd5; assign IM[19][5] = 'd20; assign IM[19][6] = 'd58; assign IM[19][7] = 'd2; 
assign IM[20][0] = 'd4; assign IM[20][1] = 'd6; assign IM[20][2] = 'd1; assign IM[20][3] = 'd48; assign IM[20][4] = 'd110; assign IM[20][5] = 'd83; assign IM[20][6] = 'd0; assign IM[20][7] = 'd89; 
assign IM[21][0] = 'd81; assign IM[21][1] = 'd40; assign IM[21][2] = 'd37; assign IM[21][3] = 'd25; assign IM[21][4] = 'd19; assign IM[21][5] = 'd52; assign IM[21][6] = 'd61; assign IM[21][7] = 'd49; 
assign IM[22][0] = 'd124; assign IM[22][1] = 'd72; assign IM[22][2] = 'd32; assign IM[22][3] = 'd41; assign IM[22][4] = 'd33; assign IM[22][5] = 'd46; assign IM[22][6] = 'd37; assign IM[22][7] = 'd69; 
assign IM[23][0] = 'd61; assign IM[23][1] = 'd117; assign IM[23][2] = 'd121; assign IM[23][3] = 'd111; assign IM[23][4] = 'd94; assign IM[23][5] = 'd18; assign IM[23][6] = 'd96; assign IM[23][7] = 'd40; 
assign IM[24][0] = 'd104; assign IM[24][1] = 'd109; assign IM[24][2] = 'd45; assign IM[24][3] = 'd87; assign IM[24][4] = 'd102; assign IM[24][5] = 'd39; assign IM[24][6] = 'd59; assign IM[24][7] = 'd90; 
assign IM[25][0] = 'd62; assign IM[25][1] = 'd76; assign IM[25][2] = 'd112; assign IM[25][3] = 'd0; assign IM[25][4] = 'd21; assign IM[25][5] = 'd70; assign IM[25][6] = 'd2; assign IM[25][7] = 'd37; 
assign IM[26][0] = 'd118; assign IM[26][1] = 'd18; assign IM[26][2] = 'd111; assign IM[26][3] = 'd56; assign IM[26][4] = 'd29; assign IM[26][5] = 'd36; assign IM[26][6] = 'd58; assign IM[26][7] = 'd122; 
assign IM[27][0] = 'd85; assign IM[27][1] = 'd87; assign IM[27][2] = 'd117; assign IM[27][3] = 'd43; assign IM[27][4] = 'd45; assign IM[27][5] = 'd103; assign IM[27][6] = 'd57; assign IM[27][7] = 'd37; 
assign IM[28][0] = 'd73; assign IM[28][1] = 'd111; assign IM[28][2] = 'd56; assign IM[28][3] = 'd95; assign IM[28][4] = 'd14; assign IM[28][5] = 'd120; assign IM[28][6] = 'd56; assign IM[28][7] = 'd45; 
assign IM[29][0] = 'd67; assign IM[29][1] = 'd93; assign IM[29][2] = 'd12; assign IM[29][3] = 'd93; assign IM[29][4] = 'd22; assign IM[29][5] = 'd0; assign IM[29][6] = 'd113; assign IM[29][7] = 'd80; 
assign IM[30][0] = 'd18; assign IM[30][1] = 'd10; assign IM[30][2] = 'd12; assign IM[30][3] = 'd116; assign IM[30][4] = 'd93; assign IM[30][5] = 'd85; assign IM[30][6] = 'd2; assign IM[30][7] = 'd2; 
assign IM[31][0] = 'd2; assign IM[31][1] = 'd36; assign IM[31][2] = 'd0; assign IM[31][3] = 'd57; assign IM[31][4] = 'd56; assign IM[31][5] = 'd121; assign IM[31][6] = 'd3; assign IM[31][7] = 'd20; 
assign IM[32][0] = 'd29; assign IM[32][1] = 'd79; assign IM[32][2] = 'd70; assign IM[32][3] = 'd10; assign IM[32][4] = 'd42; assign IM[32][5] = 'd96; assign IM[32][6] = 'd104; assign IM[32][7] = 'd49; 
assign IM[33][0] = 'd44; assign IM[33][1] = 'd10; assign IM[33][2] = 'd4; assign IM[33][3] = 'd69; assign IM[33][4] = 'd68; assign IM[33][5] = 'd10; assign IM[33][6] = 'd112; assign IM[33][7] = 'd65; 
assign IM[34][0] = 'd32; assign IM[34][1] = 'd29; assign IM[34][2] = 'd42; assign IM[34][3] = 'd48; assign IM[34][4] = 'd103; assign IM[34][5] = 'd115; assign IM[34][6] = 'd39; assign IM[34][7] = 'd12; 
assign IM[35][0] = 'd59; assign IM[35][1] = 'd113; assign IM[35][2] = 'd97; assign IM[35][3] = 'd41; assign IM[35][4] = 'd15; assign IM[35][5] = 'd127; assign IM[35][6] = 'd20; assign IM[35][7] = 'd94; 
assign IM[36][0] = 'd78; assign IM[36][1] = 'd47; assign IM[36][2] = 'd24; assign IM[36][3] = 'd111; assign IM[36][4] = 'd107; assign IM[36][5] = 'd106; assign IM[36][6] = 'd1; assign IM[36][7] = 'd6; 
assign IM[37][0] = 'd122; assign IM[37][1] = 'd68; assign IM[37][2] = 'd103; assign IM[37][3] = 'd32; assign IM[37][4] = 'd22; assign IM[37][5] = 'd27; assign IM[37][6] = 'd39; assign IM[37][7] = 'd68; 
assign IM[38][0] = 'd56; assign IM[38][1] = 'd21; assign IM[38][2] = 'd78; assign IM[38][3] = 'd7; assign IM[38][4] = 'd100; assign IM[38][5] = 'd72; assign IM[38][6] = 'd90; assign IM[38][7] = 'd71; 
assign IM[39][0] = 'd17; assign IM[39][1] = 'd3; assign IM[39][2] = 'd81; assign IM[39][3] = 'd64; assign IM[39][4] = 'd8; assign IM[39][5] = 'd13; assign IM[39][6] = 'd12; assign IM[39][7] = 'd19; 
assign IM[40][0] = 'd123; assign IM[40][1] = 'd68; assign IM[40][2] = 'd126; assign IM[40][3] = 'd76; assign IM[40][4] = 'd4; assign IM[40][5] = 'd10; assign IM[40][6] = 'd115; assign IM[40][7] = 'd119; 
assign IM[41][0] = 'd6; assign IM[41][1] = 'd22; assign IM[41][2] = 'd59; assign IM[41][3] = 'd5; assign IM[41][4] = 'd103; assign IM[41][5] = 'd86; assign IM[41][6] = 'd95; assign IM[41][7] = 'd53; 
assign IM[42][0] = 'd6; assign IM[42][1] = 'd52; assign IM[42][2] = 'd2; assign IM[42][3] = 'd61; assign IM[42][4] = 'd106; assign IM[42][5] = 'd50; assign IM[42][6] = 'd4; assign IM[42][7] = 'd79; 
assign IM[43][0] = 'd18; assign IM[43][1] = 'd94; assign IM[43][2] = 'd116; assign IM[43][3] = 'd46; assign IM[43][4] = 'd103; assign IM[43][5] = 'd1; assign IM[43][6] = 'd98; assign IM[43][7] = 'd109; 
assign IM[44][0] = 'd37; assign IM[44][1] = 'd89; assign IM[44][2] = 'd28; assign IM[44][3] = 'd2; assign IM[44][4] = 'd19; assign IM[44][5] = 'd70; assign IM[44][6] = 'd79; assign IM[44][7] = 'd121; 
assign IM[45][0] = 'd122; assign IM[45][1] = 'd66; assign IM[45][2] = 'd6; assign IM[45][3] = 'd93; assign IM[45][4] = 'd116; assign IM[45][5] = 'd7; assign IM[45][6] = 'd24; assign IM[45][7] = 'd61; 
assign IM[46][0] = 'd72; assign IM[46][1] = 'd89; assign IM[46][2] = 'd36; assign IM[46][3] = 'd100; assign IM[46][4] = 'd85; assign IM[46][5] = 'd49; assign IM[46][6] = 'd105; assign IM[46][7] = 'd29; 
assign IM[47][0] = 'd90; assign IM[47][1] = 'd125; assign IM[47][2] = 'd70; assign IM[47][3] = 'd27; assign IM[47][4] = 'd79; assign IM[47][5] = 'd95; assign IM[47][6] = 'd91; assign IM[47][7] = 'd13; 
assign IM[48][0] = 'd5; assign IM[48][1] = 'd92; assign IM[48][2] = 'd12; assign IM[48][3] = 'd121; assign IM[48][4] = 'd4; assign IM[48][5] = 'd64; assign IM[48][6] = 'd95; assign IM[48][7] = 'd56; 
assign IM[49][0] = 'd24; assign IM[49][1] = 'd34; assign IM[49][2] = 'd44; assign IM[49][3] = 'd74; assign IM[49][4] = 'd86; assign IM[49][5] = 'd110; assign IM[49][6] = 'd114; assign IM[49][7] = 'd29; 
assign IM[50][0] = 'd9; assign IM[50][1] = 'd23; assign IM[50][2] = 'd71; assign IM[50][3] = 'd49; assign IM[50][4] = 'd110; assign IM[50][5] = 'd1; assign IM[50][6] = 'd31; assign IM[50][7] = 'd21; 
assign IM[51][0] = 'd102; assign IM[51][1] = 'd37; assign IM[51][2] = 'd87; assign IM[51][3] = 'd123; assign IM[51][4] = 'd107; assign IM[51][5] = 'd104; assign IM[51][6] = 'd43; assign IM[51][7] = 'd120; 
assign IM[52][0] = 'd55; assign IM[52][1] = 'd55; assign IM[52][2] = 'd84; assign IM[52][3] = 'd66; assign IM[52][4] = 'd124; assign IM[52][5] = 'd59; assign IM[52][6] = 'd104; assign IM[52][7] = 'd91; 
assign IM[53][0] = 'd0; assign IM[53][1] = 'd73; assign IM[53][2] = 'd18; assign IM[53][3] = 'd18; assign IM[53][4] = 'd126; assign IM[53][5] = 'd9; assign IM[53][6] = 'd107; assign IM[53][7] = 'd93; 
assign IM[54][0] = 'd49; assign IM[54][1] = 'd41; assign IM[54][2] = 'd46; assign IM[54][3] = 'd42; assign IM[54][4] = 'd23; assign IM[54][5] = 'd126; assign IM[54][6] = 'd126; assign IM[54][7] = 'd77; 
assign IM[55][0] = 'd21; assign IM[55][1] = 'd121; assign IM[55][2] = 'd41; assign IM[55][3] = 'd110; assign IM[55][4] = 'd92; assign IM[55][5] = 'd35; assign IM[55][6] = 'd23; assign IM[55][7] = 'd106; 
assign IM[56][0] = 'd28; assign IM[56][1] = 'd122; assign IM[56][2] = 'd28; assign IM[56][3] = 'd74; assign IM[56][4] = 'd35; assign IM[56][5] = 'd93; assign IM[56][6] = 'd27; assign IM[56][7] = 'd116; 
assign IM[57][0] = 'd33; assign IM[57][1] = 'd56; assign IM[57][2] = 'd87; assign IM[57][3] = 'd17; assign IM[57][4] = 'd65; assign IM[57][5] = 'd112; assign IM[57][6] = 'd11; assign IM[57][7] = 'd32; 
assign IM[58][0] = 'd125; assign IM[58][1] = 'd28; assign IM[58][2] = 'd3; assign IM[58][3] = 'd42; assign IM[58][4] = 'd98; assign IM[58][5] = 'd5; assign IM[58][6] = 'd2; assign IM[58][7] = 'd17; 
assign IM[59][0] = 'd25; assign IM[59][1] = 'd55; assign IM[59][2] = 'd93; assign IM[59][3] = 'd63; assign IM[59][4] = 'd40; assign IM[59][5] = 'd89; assign IM[59][6] = 'd12; assign IM[59][7] = 'd81; 
assign IM[60][0] = 'd112; assign IM[60][1] = 'd46; assign IM[60][2] = 'd127; assign IM[60][3] = 'd104; assign IM[60][4] = 'd70; assign IM[60][5] = 'd49; assign IM[60][6] = 'd119; assign IM[60][7] = 'd114; 
assign IM[61][0] = 'd92; assign IM[61][1] = 'd115; assign IM[61][2] = 'd38; assign IM[61][3] = 'd72; assign IM[61][4] = 'd92; assign IM[61][5] = 'd123; assign IM[61][6] = 'd74; assign IM[61][7] = 'd56; 
assign IM[62][0] = 'd102; assign IM[62][1] = 'd59; assign IM[62][2] = 'd98; assign IM[62][3] = 'd96; assign IM[62][4] = 'd36; assign IM[62][5] = 'd118; assign IM[62][6] = 'd59; assign IM[62][7] = 'd84; 
assign IM[63][0] = 'd45; assign IM[63][1] = 'd119; assign IM[63][2] = 'd103; assign IM[63][3] = 'd61; assign IM[63][4] = 'd117; assign IM[63][5] = 'd33; assign IM[63][6] = 'd70; assign IM[63][7] = 'd45;


//LUT
/*integer i,i_cf,j;
always_comb begin
for (i=0;i<E/channel_fold_factor;i=i+1) begin
for (i_cf=0;i_cf<channel_fold_factor;i_cf=i_cf+1) begin
for (j=0;j<vector_fold_factor;j=j+1) begin
    //LBP_hvs[i] = IM[LBP_codes[i]]
    if ((vf_counter == j) & (cf_counter == i_cf)) begin
    case(LBP_codes[i+i_cf*E/channel_fold_factor])
        'd0 : LBP_hvs[i] = IM[0][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd1 : LBP_hvs[i] = IM[1][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd2 : LBP_hvs[i] = IM[2][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd3 : LBP_hvs[i] = IM[3][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd4 : LBP_hvs[i] = IM[4][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd5 : LBP_hvs[i] = IM[5][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd6 : LBP_hvs[i] = IM[6][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd7 : LBP_hvs[i] = IM[7][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd8 : LBP_hvs[i] = IM[8][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd9 : LBP_hvs[i] = IM[9][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd10: LBP_hvs[i] = IM[10][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd11: LBP_hvs[i] = IM[11][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd12: LBP_hvs[i] = IM[12][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd13: LBP_hvs[i] = IM[13][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd14: LBP_hvs[i] = IM[14][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd15: LBP_hvs[i] = IM[15][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd16: LBP_hvs[i] = IM[16][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd17: LBP_hvs[i] = IM[17][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd18: LBP_hvs[i] = IM[18][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd19: LBP_hvs[i] = IM[19][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd20: LBP_hvs[i] = IM[20][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd21: LBP_hvs[i] = IM[21][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd22: LBP_hvs[i] = IM[22][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd23: LBP_hvs[i] = IM[23][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd24: LBP_hvs[i] = IM[24][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd25: LBP_hvs[i] = IM[25][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd26: LBP_hvs[i] = IM[26][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd27: LBP_hvs[i] = IM[27][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd28: LBP_hvs[i] = IM[28][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd29: LBP_hvs[i] = IM[29][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd30: LBP_hvs[i] = IM[30][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd31: LBP_hvs[i] = IM[31][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd32: LBP_hvs[i] = IM[32][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd33: LBP_hvs[i] = IM[33][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd34: LBP_hvs[i] = IM[34][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd35: LBP_hvs[i] = IM[35][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd36: LBP_hvs[i] = IM[36][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd37: LBP_hvs[i] = IM[37][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd38: LBP_hvs[i] = IM[38][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd39: LBP_hvs[i] = IM[39][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd40: LBP_hvs[i] = IM[40][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd41: LBP_hvs[i] = IM[41][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd42: LBP_hvs[i] = IM[42][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd43: LBP_hvs[i] = IM[43][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd44: LBP_hvs[i] = IM[44][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd45: LBP_hvs[i] = IM[45][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd46: LBP_hvs[i] = IM[46][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd47: LBP_hvs[i] = IM[47][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd48: LBP_hvs[i] = IM[48][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd49: LBP_hvs[i] = IM[49][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd50: LBP_hvs[i] = IM[50][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd51: LBP_hvs[i] = IM[51][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd52: LBP_hvs[i] = IM[52][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd53: LBP_hvs[i] = IM[53][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd54: LBP_hvs[i] = IM[54][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd55: LBP_hvs[i] = IM[55][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd56: LBP_hvs[i] = IM[56][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd57: LBP_hvs[i] = IM[57][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd58: LBP_hvs[i] = IM[58][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd59: LBP_hvs[i] = IM[59][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd60: LBP_hvs[i] = IM[60][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd61: LBP_hvs[i] = IM[61][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd62: LBP_hvs[i] = IM[62][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	'd63: LBP_hvs[i] = IM[63][j*NB_SEGMENTS/vector_fold_factor +: NB_SEGMENTS/vector_fold_factor];
	default: LBP_hvs[i] = 'd0;
    endcase
  end
end
end
end
end*/

always_comb begin
  for (int i = 0; i < E/channel_fold_factor; i = i + 1) begin
    
    // 1. Calculate dynamic indices
    automatic int code_idx   = i + (cf_counter * (E / channel_fold_factor));
    automatic int slice_base = vf_counter * (NB_SEGMENTS / vector_fold_factor);
    
    // 2. Fetch the target index from LBP_codes
    automatic int target_im_idx = LBP_codes[code_idx];

    // 3. Assign the value with bounds checking
    if (target_im_idx <= 63) begin
      // The width portion MUST be a constant expression directly in the brackets
      LBP_hvs[i] = IM[target_im_idx][slice_base +: (NB_SEGMENTS / vector_fold_factor)];
    end else begin
      LBP_hvs[i] = 'd0; 
    end

  end
end


endmodule

module BINDING_o1 # (
    parameter integer D = 1024,
    //parameter real p_sparse = 0.0078125,
    parameter integer E = 64,
    parameter integer NB_SEGMENTS = 8,
    parameter integer LENGTH_SEGMENT = 128,
    parameter integer LBP_LENGTH = 6,
    parameter integer vector_fold_factor = 1, //either 1,2,4...NB_SEGMENTS
    parameter integer channel_fold_factor = 1 //either 1,2,4...E (if E is a power of 2)
)
(
    input  [E/channel_fold_factor-1:0][(NB_SEGMENTS/vector_fold_factor)-1:0][$clog2(LENGTH_SEGMENT)-1:0] LBP_hvs,
    input  [$clog2(vector_fold_factor)-1:0] vf_counter,
    input  [$clog2(channel_fold_factor)-1:0] cf_counter,
    output [E/channel_fold_factor-1:0][0:LENGTH_SEGMENT*NB_SEGMENTS/vector_fold_factor-1] bind_outputs

);
logic [$clog2(LENGTH_SEGMENT):0][E/channel_fold_factor-1:0][0:LENGTH_SEGMENT*NB_SEGMENTS/vector_fold_factor-1] intermediate_bind_outputs;
assign bind_outputs = intermediate_bind_outputs[$clog2(LENGTH_SEGMENT)];


logic [E-1:0][0:LENGTH_SEGMENT*NB_SEGMENTS-1] EM;

assign EM[0] = 1024'b0000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000;
assign EM[1] = 1024'b0000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000;
assign EM[2] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000;
assign EM[3] = 1024'b0000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[4] = 1024'b0000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[5] = 1024'b0000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000;
assign EM[6] = 1024'b0000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000;
assign EM[7] = 1024'b0000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000;
assign EM[8] = 1024'b0001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000101000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000;
assign EM[9] = 1024'b0000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000;
assign EM[10] = 1024'b0000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[11] = 1024'b0000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[12] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100;
assign EM[13] = 1024'b0000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[14] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000;
assign EM[15] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000;
assign EM[16] = 1024'b0000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000;
assign EM[17] = 1024'b0100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[18] = 1024'b0000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[19] = 1024'b0000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000;
assign EM[20] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000;
assign EM[21] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[22] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000;
assign EM[23] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[24] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000;
assign EM[25] = 1024'b0001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[26] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[27] = 1024'b0000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000;
assign EM[28] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[29] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[30] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000;
assign EM[31] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000;
assign EM[32] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[33] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000;
assign EM[34] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[35] = 1024'b0000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[36] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[37] = 1024'b0000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000;
assign EM[38] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000;
assign EM[39] = 1024'b0000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000;
assign EM[40] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[41] = 1024'b0000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000;
assign EM[42] = 1024'b0000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000;
assign EM[43] = 1024'b0000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[44] = 1024'b0000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[45] = 1024'b0000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[46] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000;
assign EM[47] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000;
assign EM[48] = 1024'b0000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000;
assign EM[49] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000;
assign EM[50] = 1024'b0000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000;
assign EM[51] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000;
assign EM[52] = 1024'b0100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[53] = 1024'b0000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[54] = 1024'b0000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[55] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[56] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000;
assign EM[57] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[58] = 1024'b0000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[59] = 1024'b0000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000;
assign EM[60] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000;
assign EM[61] = 1024'b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000;
assign EM[62] = 1024'b0000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000;
assign EM[63] = 1024'b0000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;

//set EM as inputs to the binding
integer k,k50,k_cf;
always_comb begin
for (k=0;k<vector_fold_factor;k=k+1) begin
for (k_cf=0;k_cf<channel_fold_factor;k_cf=k_cf+1) begin
    if ((vf_counter == k) & (cf_counter == k_cf)) begin
      for (k50=0;k50<E/channel_fold_factor;k50=k50+1) begin
          intermediate_bind_outputs[0][k50] = EM[k50+k_cf*E/channel_fold_factor][k*NB_SEGMENTS/vector_fold_factor*LENGTH_SEGMENT +: NB_SEGMENTS/vector_fold_factor*LENGTH_SEGMENT];
      end
    end
end
end
end

//All the E/channel_fold_factor*NB_SEGMENTS/vector_fold_factor binding units
genvar i2,i3,i4;
generate 
    for (i4=0;i4<E/channel_fold_factor;i4=i4+1) begin: binding_i4_gen
    for (i3=0;i3<NB_SEGMENTS/vector_fold_factor;i3=i3+1) begin: binding_i3_gen
    for (i2=0;i2<$clog2(LENGTH_SEGMENT);i2=i2+1) begin: binding_i2_gen
	select_if_shift_o1 #(.LENGTH_SEGMENT(LENGTH_SEGMENT),.SHIFT(2**(i2))) select_if_shift0 (.segment(intermediate_bind_outputs[i2][i4][LENGTH_SEGMENT*i3:LENGTH_SEGMENT*(i3+1)-1]), .enable(LBP_hvs[i4][i3][i2]), .result(intermediate_bind_outputs[i2+1][i4][LENGTH_SEGMENT*i3:LENGTH_SEGMENT*(i3+1)-1]));
    end
    end
    end
endgenerate

endmodule

module BUNDLING_SPACE_o1 # (
    parameter integer D = 1024,
    //parameter real p_sparse = 0.0078125,
    parameter integer E = 64,
    parameter integer NB_SEGMENTS = 8,
    parameter integer LENGTH_SEGMENT = 128,
    parameter integer LBP_LENGTH = 6,
    parameter integer vector_fold_factor = 1, //either 1,2,4...NB_SEGMENTS
    parameter integer channel_fold_factor = 1 //either 1,2,4...E (if E is a power of 2)
)
(
    input  clk_i,
    input  arst_n_in,
    input  [E/channel_fold_factor-1:0][0:LENGTH_SEGMENT*NB_SEGMENTS/vector_fold_factor-1] bind_outputs,
    input  [$clog2(channel_fold_factor)-1:0] cf_counter,
    output logic bundled_hv [0:LENGTH_SEGMENT*NB_SEGMENTS/vector_fold_factor-1]
);

logic [E/channel_fold_factor-1:0] bound_hvs [0:LENGTH_SEGMENT*NB_SEGMENTS/vector_fold_factor-1];
logic bundled_hv_e [0:LENGTH_SEGMENT*NB_SEGMENTS/vector_fold_factor-1];
logic extra_ch_fold_reg [0:LENGTH_SEGMENT*NB_SEGMENTS/vector_fold_factor-1];

//By switching the packed and unpacked parts, can have very simple OR-tree module
genvar i10,i11;
generate
for (i11=0;i11<D/vector_fold_factor;i11=i11+1) begin
for (i10=0;i10<E/channel_fold_factor;i10=i10+1) begin
    assign bound_hvs[i11][i10] = bind_outputs[i10][i11];
end
end
endgenerate


//Bundle in space
OR_tree_o1 #(.D(D/vector_fold_factor),.NB_INPUTS(E/channel_fold_factor)) OR_tree0 (.input_hvs(bound_hvs), .output_hv(bundled_hv_e));

integer i_bun;
always_comb begin
for (i_bun=0;i_bun<LENGTH_SEGMENT*NB_SEGMENTS/vector_fold_factor;i_bun=i_bun+1) begin
    bundled_hv[i_bun] = bundled_hv_e[i_bun];
end
end

integer i_bun2,i_bun3;
always @(posedge clk_i) begin
    for (i_bun2=0;i_bun2<LENGTH_SEGMENT*NB_SEGMENTS/vector_fold_factor;i_bun2=i_bun2+1) begin
    for (i_bun3=0;i_bun3<channel_fold_factor;i_bun3=i_bun3+1) begin
    if (arst_n_in == 0) begin
	extra_ch_fold_reg[i_bun2] <= 'd0;
    end
    else if (i_bun3 == cf_counter)
        extra_ch_fold_reg[i_bun2] <= (i_bun3 == (channel_fold_factor-1)) ? 'd0:(bundled_hv_e[i_bun2] | extra_ch_fold_reg[i_bun2]);
    end
    end
end


endmodule


module BUNDLING_TIME_o1 # (
    parameter integer D = 1024,
    //parameter real p_sparse = 0.0078125,
    parameter integer E = 64,
    parameter integer NB_TO_BUNDLE_IN_TIME = 256,
    parameter integer NB_SEGMENTS = 8,
    parameter integer LENGTH_SEGMENT = 128,
    parameter integer LBP_LENGTH = 6,
    parameter integer vector_fold_factor = 1, //either 1,2,4...NB_SEGMENTS
    parameter integer channel_fold_factor = 1 //either 1,2,4...E (if E is a power of 2)
)
(
    input  clk_i,
    input  arst_n_in,
    input  bundled_hv [0:LENGTH_SEGMENT*NB_SEGMENTS/vector_fold_factor-1],
    input  [$clog2(vector_fold_factor)-1:0] vf_counter,
    input  [$clog2(channel_fold_factor)-1:0] cf_counter,
    output logic [0:D-1] result_hv,
    output logic [$clog2(NB_TO_BUNDLE_IN_TIME):0] counter_bund_time,
    output logic [$clog2(NB_TO_BUNDLE_IN_TIME)-1:0] counter_array [0:D-1]
);




integer i15,j2,j3,j4;
always @(posedge clk_i or negedge arst_n_in) begin
	if (arst_n_in == 0) begin
	    counter_bund_time <= 9'd0;
    	    for (i15=0;i15<D/vector_fold_factor;i15=i15+1) begin
	    counter_array[i15] <= 8'd0;
            for (j4=1;j4<vector_fold_factor;j4=j4+1) begin
		counter_array[i15+j4*D/vector_fold_factor] <= 8'd0;
	    end
	    end
	end
	else begin
   	for (i15=0;i15<D/vector_fold_factor;i15=i15+1) begin
	   if ((counter_bund_time == NB_TO_BUNDLE_IN_TIME) & (cf_counter == channel_fold_factor-1)) begin
            for (j3=0;j3<vector_fold_factor;j3=j3+1) begin
		if (vf_counter == j3) begin
		    counter_array[i15+j3*D/vector_fold_factor] <= bundled_hv[i15];
		    if (vf_counter == (vector_fold_factor-1)) begin
			counter_bund_time <= 1;
		    end
		end
	    end
	end
	else if (cf_counter == channel_fold_factor-1) begin
    	    for (j2=0;j2<vector_fold_factor;j2=j2+1) begin
	        if (vf_counter == j2) begin
                    counter_array[i15+j2*D/vector_fold_factor] <= (bundled_hv[i15]==1) ? counter_array[i15+j2*D/vector_fold_factor]+1:counter_array[i15+j2*D/vector_fold_factor];
	        end
	    end
	    counter_bund_time <= (vf_counter == (vector_fold_factor-1)) ? counter_bund_time+1:counter_bund_time;
	end
    end
    end
end


genvar i20;
generate
for (i20=0;i20<D;i20=i20+1) begin
always @(posedge clk_i) begin
    if ((counter_bund_time == NB_TO_BUNDLE_IN_TIME) & (vf_counter == 0) & (cf_counter == channel_fold_factor-1)) begin
	result_hv[i20] <= (counter_array[i20] > 8'h82) ? (1'b1):(1'b0); 
    end
end
end
endgenerate

endmodule


module SIMILARITY_SEARCH_o1 # (
    parameter integer D = 1024,
    //parameter real p_sparse = 0.0078125,
    parameter integer E = 64,
    parameter integer NB_TO_BUNDLE_IN_TIME = 256,
    parameter integer NB_SEGMENTS = 8,
    parameter integer LENGTH_SEGMENT = 128,
    parameter integer LBP_LENGTH = 6,
    parameter integer NB_CLASSES = 2,
    parameter integer vector_fold_factor = 1, //either 1,2,4...NB_SEGMENTS
    parameter integer channel_fold_factor = 1 //either 1,2,4...E (if E is a power of 2)
)
(
    input  clk_i,
    input  arst_n_in,
    input  [0:D-1] result_hv,
    input  [0:D-1] ictal_hv,
    input  [0:D-1] interictal_hv,
    input  [$clog2(vector_fold_factor)-1:0] vf_counter,
    input  [$clog2(channel_fold_factor)-1:0] cf_counter,
    input  [$clog2(NB_TO_BUNDLE_IN_TIME):0] counter_bund_time,
    output classification,
    output classification_ready,
    output [$clog2(D):0] ictal_sim,
    output [$clog2(D):0] interictal_sim
);

logic [0:D-1] sim_vector;
logic [$clog2(D):0] sim_score;
logic [$clog2(D):0] sim_score_saved;
assign interictal_sim = sim_score_saved;
assign ictal_sim = sim_score;
logic [$clog2(NB_CLASSES):0] sim_counter;


always_comb begin
    if (sim_counter == 0) begin
	sim_vector = result_hv & interictal_hv;
    end
    else if (sim_counter == 1) begin
	sim_vector = result_hv & ictal_hv;
    end
end
adder_tree_o1 # ( .NB_INPUTS(D) ) adder_tree0 ( .input_data(sim_vector), .output_data(sim_score) );

always @(posedge clk_i) begin
    if ((counter_bund_time == NB_TO_BUNDLE_IN_TIME) & (vf_counter == 0) & (cf_counter == channel_fold_factor-1)) begin
	sim_counter <= 0;
    end
    else if (sim_counter == 0) begin
	sim_score_saved <= sim_score;
	sim_counter <= 1;
    end
    else if (sim_counter == 1) begin
	sim_counter <= 2;
    end
end
assign classification = (sim_score > sim_score_saved) ? 1 : 0;
assign classification_ready = (sim_counter == 1) ? 1 : 0;

endmodule


module select_if_shift_o1 # (
    parameter LENGTH_SEGMENT = 128,
    parameter SHIFT = 1
)
(
    input  [LENGTH_SEGMENT-1:0] segment,
    input  enable,
    output [LENGTH_SEGMENT-1:0] result
);
logic [LENGTH_SEGMENT-1:0] temp;
assign result = (enable) ? temp : segment;

always_comb begin
    temp[LENGTH_SEGMENT-1:SHIFT] = segment[LENGTH_SEGMENT-1-SHIFT:0];
    temp[SHIFT-1:0] = segment[LENGTH_SEGMENT-1:LENGTH_SEGMENT-SHIFT];
end

endmodule


module OR_tree_o1 # (
    parameter D = 1024,
    parameter NB_INPUTS = 64
)
(
    input [NB_INPUTS-1:0] input_hvs [D-1:0],
    output output_hv [D-1:0]
);

genvar i;
generate 
for (i=0;i<D;i=i+1) begin
    assign output_hv[i] = |input_hvs[i];
end
endgenerate 

endmodule

module adder_tree_o1 # (
    parameter NB_INPUTS,
    parameter NB_STAGES = $clog2(NB_INPUTS)
)
(
    input logic [NB_INPUTS-1:0] input_data,
    output logic [NB_STAGES:0] output_data
);

logic [NB_STAGES-1:0][(NB_INPUTS/2)-1:0][NB_STAGES:0] data; //3th part is output_data_width

genvar stage_nb, adder_nb;
generate
for (stage_nb=0;stage_nb<NB_STAGES;stage_nb=stage_nb+1) begin: stage_gen
    localparam nb_outputs_stage = NB_INPUTS >> (stage_nb+1);
    localparam data_width_stage = stage_nb+2;

    if (stage_nb == 0) begin
	for (adder_nb=0;adder_nb<nb_outputs_stage;adder_nb=adder_nb+1) begin: first_adder_gen
	    always_comb begin
		data[stage_nb][adder_nb][data_width_stage-1:0] = input_data[adder_nb*2] + input_data[adder_nb*2+1];
	    end
	end
    end
    else begin
	for (adder_nb=0;adder_nb<nb_outputs_stage;adder_nb=adder_nb+1) begin: adder_gen
	    always_comb begin
		data[stage_nb][adder_nb][data_width_stage-1:0] = 
			data[stage_nb-1][adder_nb*2][(data_width_stage-1)-1:0] + 
			data[stage_nb-1][adder_nb*2+1][(data_width_stage-1)-1:0];
	    end
        end
    end
end
endgenerate

assign output_data = data[NB_STAGES-1][0];

endmodule
