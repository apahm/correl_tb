`timescale 1 ns/1 ps

module ham_tb();

reg aclk_reg;
reg aresetn_reg;

wire           FFT_CONFIG_tready;
wire           IFFT_CONFIG_tready;
wire           S_AXIS_IM_tready;
wire           S_AXIS_RE_tready;
wire           S_AXIS_SF_tready;

wire           SF_IM_tready;
wire           SF_RE_tready;

reg [15:0]     FFT_CONFIG_tdata_reg;
reg            FFT_CONFIG_tvalid_reg;

reg [23:0]     IFFT_CONFIG_tdata_reg;
reg            IFFT_CONFIG_tvalid_reg;

reg [15:0]     S_AXIS_IM_tdata_reg;
reg            S_AXIS_IM_tlast_reg;
reg            S_AXIS_IM_tvalid_reg;

reg [15:0]     S_AXIS_RE_tdata_reg;
reg            S_AXIS_RE_tlast_reg;
reg            S_AXIS_RE_tvalid_reg;

reg [63:0]     S_AXIS_SF_tdata_reg;
reg            S_AXIS_SF_tvalid_reg;

reg [15:0]     SF_RE_tdata_reg;
reg [15:0]     SF_IM_tdata_reg;

reg            SF_RE_tvalid_reg;
reg            SF_IM_tvalid_reg;


design_1_wrapper
design_1_wrapper_inst
(
   // Настройка БПФ ядра  
   .FFT_CONFIG_tdata(FFT_CONFIG_tdata_reg),
   .FFT_CONFIG_tready(FFT_CONFIG_tready),
   .FFT_CONFIG_tvalid(FFT_CONFIG_tvalid_reg),

   // Настройка ОБПФ ядра
   .IFFT_CONFIG_tdata(IFFT_CONFIG_tdata_reg),
   .IFFT_CONFIG_tready(IFFT_CONFIG_tready),
   .IFFT_CONFIG_tvalid(IFFT_CONFIG_tvalid_reg),

   .M_AXIS_IM_tdata(),
   .M_AXIS_IM_tready(1'b1),
   .M_AXIS_IM_tvalid(),
   
   .M_AXIS_RE_tdata(),
   .M_AXIS_RE_tlast(),
   .M_AXIS_RE_tready(1'b1),
   .M_AXIS_RE_tvalid(),
   
   .S_AXIS_IM_tdata(S_AXIS_IM_tdata_reg),
   .S_AXIS_IM_tlast(S_AXIS_IM_tlast_reg),
   .S_AXIS_IM_tvalid(S_AXIS_RE_tvalid_reg),
   .S_AXIS_IM_tready(S_AXIS_IM_tready),
   
   .S_AXIS_RE_tdata(S_AXIS_RE_tdata_reg),
   .S_AXIS_RE_tlast(S_AXIS_RE_tlast_reg),
   .S_AXIS_RE_tvalid(S_AXIS_RE_tvalid_reg),
   .S_AXIS_RE_tready(S_AXIS_RE_tready),
    
   .SF_IM_tdata(SF_IM_tdata_reg),
   .SF_IM_tready(SF_IM_tready),
   .SF_IM_tvalid(SF_IM_tvalid_reg),
   
   .SF_RE_tdata(SF_RE_tdata_reg),
   .SF_RE_tready(SF_RE_tready),
   .SF_RE_tvalid(SF_RE_tvalid_reg),

   .aclk(aclk_reg),
   .aresetn(aresetn_reg)
);


initial begin
   aclk_reg = 0;
   forever begin
      #5;
      aclk_reg = !aclk_reg;
   end
end

integer fd_recv_re;
integer fd_recv_im;  

integer fd_recv_int16_re;
integer fd_recv_int16_im;  

reg [15:0] lfm_reg;

reg [15:0] lfm_re_shift_reg [1023:0];
reg [15:0] lfm_im_shift_reg [1023:0];

reg [15:0] lfm_re_int16_reg [1023:0];
reg [15:0] lfm_im_int16_reg [1023:0];

initial begin
   $display("Read LFM from files ...");

   fd_recv_re = $fopen("lfm_int16_shift_re.dat","rb");
   fd_recv_im = $fopen("lfm_int16_shift_im.dat","rb");
   $fread(lfm_re_shift_reg, fd_recv_re);
   $fread(lfm_im_shift_reg, fd_recv_im);

   fd_recv_int16_re = $fopen("lfm_re_int16.dat","rb");
   fd_recv_int16_im = $fopen("lfm_im_int16.dat","rb");
   $fread(lfm_re_int16_reg, fd_recv_int16_re);
   $fread(lfm_im_int16_reg, fd_recv_int16_im);

   $display("Begin verification of correl function windows ...");
   aresetn_reg = 0;
   S_AXIS_SF_tdata_reg = 64'b0;
   S_AXIS_SF_tvalid_reg = 1'b0;

   S_AXIS_IM_tdata_reg = 16'b0;
   S_AXIS_IM_tlast_reg = 1'b0;
   S_AXIS_IM_tvalid_reg = 1'b0;

   S_AXIS_RE_tdata_reg = 16'b0;
   S_AXIS_RE_tlast_reg = 1'b0;
   S_AXIS_RE_tvalid_reg = 1'b0;

   FFT_CONFIG_tdata_reg = 16'b0;
   FFT_CONFIG_tvalid_reg = 1'b0;
   IFFT_CONFIG_tdata_reg = 24'b0;
   IFFT_CONFIG_tvalid_reg = 1'b0;

   SF_IM_tvalid_reg = 0;
   SF_RE_tvalid_reg = 0;
   SF_RE_tdata_reg = 16'b0;
   SF_IM_tdata_reg = 16'b0;
   #100;
   aresetn_reg = 1;
   #250;
   
   // Задаём настройки для БПФ
   $display("FFT_CONFIG_tready event");
   wait(FFT_CONFIG_tready);
   FFT_CONFIG_tdata_reg = {15'b0, 1'b1};
   @(posedge aclk_reg);
   FFT_CONFIG_tvalid_reg = 1'b1;
   @(posedge aclk_reg);
   FFT_CONFIG_tvalid_reg = 1'b0;

   // Задаём настройки для ОБПФ
   $display("IFFT_CONFIG_tdata_reg event");
   wait(IFFT_CONFIG_tready);
   IFFT_CONFIG_tdata_reg = {18'b0 , 1'b0, 5'b01010};
   @(posedge aclk_reg);
   IFFT_CONFIG_tvalid_reg = 1'b1;
   @(posedge aclk_reg);
   IFFT_CONFIG_tvalid_reg = 1'b0;

   // Задаём опорную функцию для сжатия
   // Задаём принятый сигнал для сжатия

   wait(SF_RE_tready);
   wait(SF_IM_tready);
   wait(S_AXIS_RE_tready);
   @(posedge aclk_reg);
   #9.9;
   SF_RE_tdata_reg = lfm_re_int16_reg[0];
   SF_IM_tdata_reg = lfm_im_int16_reg[0];
   S_AXIS_RE_tdata_reg = lfm_re_shift_reg[0];
   S_AXIS_IM_tdata_reg = lfm_im_shift_reg[0]; 
   SF_IM_tvalid_reg = 1;
   SF_RE_tvalid_reg = 1;
   S_AXIS_RE_tvalid_reg = 1;
   @(posedge aclk_reg);
   for (integer i = 1; i < 1024; i = i + 1) begin
      @(posedge aclk_reg);
      SF_RE_tdata_reg = lfm_re_int16_reg[i];
      SF_IM_tdata_reg = lfm_im_int16_reg[i];
      S_AXIS_RE_tdata_reg = lfm_re_shift_reg[i];
      S_AXIS_IM_tdata_reg = lfm_im_shift_reg[i]; 
   end
   #9.9;   
   SF_RE_tvalid_reg = 0;
   SF_IM_tvalid_reg = 0;
   S_AXIS_RE_tvalid_reg = 0;

   $display("TEST PASSED");
end

endmodule