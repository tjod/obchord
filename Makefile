#INCS = -I/usr/include/postgresql/9.1/server -I/usr/include/pgsql/server
#INCS = -I/usr/include/postgresql/8.4/server -I/usr/include/pgsql/server
INCS = -I`pg_config --includedir-server`
PG_LIB = `pg_config --pkglibdir`
#CC = /usr/bin/gcc -fpic

varbit.so: varbit.c
	$(CC) $(INCS) -fPIC -shared -o $@ $?

install:
	cp varbit.so ${PG_LIB}
