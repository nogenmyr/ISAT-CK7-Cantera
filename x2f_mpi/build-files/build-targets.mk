# Author: Varun Hiremath <vh63@cornell.edu>
# Date: Tue, 29 Nov 2011 20:37:48 -0600

$(FILES):
	$(LN) $(SOURCE_DIR)/$@ .

%.o: %.f90
	$(FC) $(FCFLAGS) $(IFLAGS) -c $<

%.o: %.f
	$(FC) $(FCFLAGS) $(IFLAGS) -c $<

%.o: %.c
	$(CC) $(CCFLAGS) $(IFLAGS) -c $<

$(LIB_NAME).so: $(OBJS)
	$(FC) -shared -o $@ $(OBJS) $(LDFLAGS)

$(LIB_NAME).a: $(OBJS)
	$(AR) $@ $(OBJS)

$(NAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LINK_LIBS) $(LDFLAGS) $(LINK_SYS_LIBS)
ifneq ($(BUILD_TYPE), debug)
	strip $(NAME)
endif

debug: clean
	make BUILD_TYPE=debug

clean::
	$(RM) $(CREATE_LIBS) $(NAME) $(OBJS) *.mod
